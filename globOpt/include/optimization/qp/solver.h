#ifndef __GF2_SOLVER_H__
#define __GF2_SOLVER_H__

//////////////
/// Solver
//////////////

#ifdef GF2_USE_GUROBI
#   include "optimization/qp/gurobiOpt.h"
#endif
#include "globfit2/primitives/pointPrimitive.h"
#include "globfit2/primitives/linePrimitive2.h" // remove, if typedef is moved
#include "qcqpcpp/io/io.h"                      // read/writeSparseMatrix

namespace GF2 {

struct SolverParams
{
        int n_points = 50;
}; // ... struct SolverParams

class Solver
{
    public:
        typedef float                                       Scalar;
        typedef Eigen::Matrix<Scalar,3,1>                   Vector;
        typedef LinePrimitive2                              PrimitiveT;
        typedef PointPrimitive                              PointPrimitiveT;
        typedef std::vector<std::vector<PrimitiveT> >       PrimitiveContainerT;
        typedef std::vector<PointPrimitiveT>                PointContainerT;
        typedef Eigen::SparseMatrix<Scalar,Eigen::RowMajor> SparseMatrix;

        //static inline int show       ( int argc, char** argv );
#if WITH_SAMPLE_INPUT
        static inline int sampleInput( int argc, char** argv );
#endif // WITH_SAMPLE_INPUT
        static inline int generateCli   ( int argc, char** argv );
        //static inline int formulate  ( int argc, char** argv );
        static inline int solve      ( int argc, char** argv );
        static inline int datafit    ( int argc, char** argv );

        //static inline int run        ( std::string img_path, Scalar const scale, std::vector<Scalar> const& angles, int argc, char** argv ) __attribute__ ((deprecated));

        static inline Eigen::Matrix<Scalar,3,1> checkSolution( std::vector<Scalar>       const& x
                                                             , SparseMatrix              const& qo
                                                             , SparseMatrix              const& Qo
                                                             , SparseMatrix              const& A
                                                             , Eigen::Matrix<Scalar,3,1> const& weights );
}; // ... cls Solver
} // ... ns gf2


//__________________________________HPP__________________________________________________

#include "Eigen/Sparse"

#ifdef GF2_USE_PCL
//#   include <pcl/visualization/pcl_visualizer.h>
#   include <pcl/console/parse.h>
#endif // GF2_USE_PCL

#ifdef GF2_WITH_MOSEK
#   include "qcqpcpp/mosekOptProblem.h"
#endif // GF2_WITH_MOSEK
#ifdef GF2_WITH_BONMIN
#   include "qcqpcpp/bonminOptProblem.h"
#endif
#ifdef GF2_WITH_GUROBI
//#   include "qcqpcpp/gurobiOptProblem.h"
#endif

#include "globfit2/util/diskUtil.hpp"
#include "globfit2/util/util.hpp"                     // timestamp2Str

#include "globfit2/primitives/pointPrimitive.h"
#include "globfit2/primitives/linePrimitive2.h"

#include "globfit2/io/io.h"
#include "globfit2/ground_truth/gtCreator.h"
#include "globfit2/optimization/candidateGenerator.h"
#include "globfit2/optimization/energyFunctors.h"     // PointLineDistanceFunctor,
#include "globfit2/optimization/problemSetup.h"       // everyPatchNeedsDirection()

namespace GF2
{

#if GF2_WITH_SAMPLE_INPUT
/**! \brief      Step 0. Takes an image, and samples it to a pointcloud. Saves points to img_path.parent_path/cloud.ply.
*    \param argc Number of CLI arguments.
*    \param argv Vector of CLI arguments.
*    \return     EXIT_SUCCESS.
*/
int
Solver::sampleInput( int argc, char** argv )
{
    std::string img_path;
    bool valid_input = true;
    if ( (pcl::console::parse_argument(argc, argv, "--img"       , img_path  ) < 0) )
    {
        std::cerr << "[" << __func__ << "]: " << "--img is compulsory" << std::endl;
        valid_input = false;
    }

    int         n_points    = 200;
    float       scene_size  = 1.f,
                noise       = 0.015f,
                filter_size = 0.0075f;
    std::string out_dir = ".";
    pcl::console::parse_argument(argc, argv, "-N"           , n_points   );
    pcl::console::parse_argument(argc, argv, "--out_dir"    , out_dir    );
    pcl::console::parse_argument(argc, argv, "--scene-size" , scene_size );
    pcl::console::parse_argument(argc, argv, "--noise"      , noise      );
    pcl::console::parse_argument(argc, argv, "--filter-cell", filter_size);
    {
         std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt --sample-input"
                   << " --img "          << img_path
                   << " [-N "            << n_points << "]"
                   << " [--out_dir "     << out_dir << "]"
                   << " [--scene-size "  << scene_size << "]"
                   << " [--noise "       << noise << "]"
                   << " [--filter-cell " << filter_size << "]"
                   << std::endl;

         if ( !valid_input )
             return EXIT_FAILURE;
    }

    // sample image
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud( new pcl::PointCloud<pcl::PointXYZRGB>() );
    GTCreator::sampleImage( cloud, img_path, n_points, noise, scene_size, filter_size, /* sensor_origin: */ NULL );

    // create output directory
    std::string img_name = boost::filesystem::path( img_path ).filename().stem().string();
    std::stringstream ss;
    ss << out_dir << "/" << img_name << "_noise_" << noise;
    std::string out_path = ss.str();
    if ( !boost::filesystem::exists(out_path) )
        boost::filesystem::create_directory( out_path );

    std::string cloud_path = out_path + "/cloud.ply";
    util::saveBackup( cloud_path );
    pcl::io::savePLYFileASCII( cloud_path, *cloud );
    std::cout << "[" << __func__ << "]: " << "saved " << cloud_path << std::endl;

    return EXIT_SUCCESS;
} // ...Solver::sampleInput()
#endif // GF2_WITH_SAMPLE_INPUT

//! \brief                  Step 1. Generates primitives from a cloud. Reads "cloud.ply" and saves "candidates.csv".
//! \param argc             Contains --cloud cloud.ply, and --scale scale.
//! \param argv             Contains --cloud cloud.ply, and --scale scale.
//! \return                 EXIT_SUCCESS.
int
Solver::generateCli( int    argc
                   , char** argv )
{
    typedef typename PointContainerT::value_type PointPrimitiveT;
    int err = EXIT_SUCCESS;

    CandidateGeneratorParams<Scalar> generatorParams;
    std::string                 cloud_path              = "./cloud.ply";
    std::vector<Scalar>         angle_gens              = { Scalar(90.) };
    std::string                 mode_string             = "representative_sqr";
    std::vector<std::string>    mode_opts               = { "representative_sqr"
#if GF2_WITH_FULL_LINKAGE
                                                            "full_min", "full_max", "squared_min", "representative_min",
#endif // GF2_WITH_FULL_LINKAGE
                                                          };
    //std::string                 patch_refit_mode_string = "avg_dir";
    //std::vector<std::string>    patch_refit_mode_opts   = { "spatial", "avg_dir" };
    std::string                 input_prims_path         = "patches.csv";
    std::string                 associations_path       = "points_primitives.csv";

    // parse input
    if ( err == EXIT_SUCCESS )
    {
        bool valid_input = true;

        // cloud
        if ( (pcl::console::parse_argument( argc, argv, "--cloud", cloud_path) < 0)
             && !boost::filesystem::exists( cloud_path ) )
        {
            std::cerr << "[" << __func__ << "]: " << "--cloud does not exist: " << cloud_path << std::endl;
            valid_input = false;
        }

        // scale
        if ( (pcl::console::parse_argument( argc, argv, "--scale", generatorParams.scale) < 0) && (pcl::console::parse_argument( argc, argv, "-sc", generatorParams.scale) < 0) )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale is compulsory" << std::endl;
            valid_input = false;
        }

        if (    (pcl::console::parse_argument( argc, argv, "-p", input_prims_path) < 0)
             && (pcl::console::parse_argument( argc, argv, "--prims", input_prims_path) < 0)
             && (!boost::filesystem::exists(input_prims_path)) )
        {
            std::cerr << "[" << __func__ << "]: " << "-p or --prims is compulsory" << std::endl;
            valid_input = false;
        }

        if (    (pcl::console::parse_argument( argc, argv, "-a", associations_path) < 0)
             && (pcl::console::parse_argument( argc, argv, "--assoc", associations_path) < 0)
             && (!boost::filesystem::exists(associations_path)) )
        {
            std::cerr << "[" << __func__ << "]: " << "-a or --assoc is compulsory" << std::endl;
            valid_input = false;
        }

        pcl::console::parse_argument( argc, argv, "--angle-limit", generatorParams.angle_limit );
        pcl::console::parse_argument( argc, argv, "-al", generatorParams.angle_limit );
        pcl::console::parse_argument( argc, argv, "--angle-limit-div", generatorParams.angle_limit_div );
        pcl::console::parse_argument( argc, argv, "-ald", generatorParams.angle_limit_div );
        pcl::console::parse_argument( argc, argv, "--patch-dist-limit", generatorParams.patch_dist_limit_mult ); // gets multiplied by scale
        pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens );
        pcl::console::parse_argument( argc, argv, "--patch-pop-limit", generatorParams.patch_population_limit );

        // patchDistMode
        pcl::console::parse_argument( argc, argv, "--mode", mode_string );
        generatorParams.parsePatchDistMode( mode_string );
        // refit
        if ( pcl::console::find_switch( argc, argv, "--patch-refit" ) )
        {
            std::cerr << "[" << __func__ << "]: " << "--patch-refit option has been DEPRECATED. exiting." << std::endl;
            return EXIT_FAILURE;
        }
        //pcl::console::parse_argument( argc, argv, "--patch-refit", patch_refit_mode_string );
        //generatorParams.parseRefitMode( patch_refit_mode_string );

        // small_mode
        {
            int small_mode = 0;
            pcl::console::parse_argument( argc, argv, "--small-mode", small_mode );
            generatorParams.small_mode = static_cast<CandidateGeneratorParams<Scalar>::SmallPatchesMode>( small_mode );
        }

        // print usage
        {
            std::cerr << "[" << __func__ << "]: " << "Usage:\t " << argv[0] << " --generate \n";
            std::cerr << "\t --cloud " << cloud_path << "\n";
            std::cerr << "\t -sc,--scale " << generatorParams.scale << "\n";
            std::cerr << "\t -p,--prims" << input_prims_path << "\n";
            std::cerr << "\t -a,--assoc" << associations_path << "\n";

            // linkage mode (full_min, full_max, squared_min, repr_min)
            std::cerr << "\t [--mode *" << generatorParams.printPatchDistMode() << "*\t";
            for ( size_t m = 0; m != mode_opts.size(); ++m )
                std::cerr << "|" << mode_opts[m];
            std::cerr << "]\n";

            // patch refit mode (spatial, avg_dir)
//            std::cerr << "\t [--patch-refit *" << generatorParams.printRefitMode() << "*\t";
//            for ( size_t m = 0; m != patch_refit_mode_opts.size(); ++m )
//                std::cerr << "|" << patch_refit_mode_opts[m];
//            std::cerr << "]\n";

            std::cerr << "\t [-al,--angle-limit " << generatorParams.angle_limit << "]\n";
            std::cerr << "\t [-ald,--angle-limit-div " << generatorParams.angle_limit_div << "]\n";
            std::cerr << "\t [--patch-dist-limit " << generatorParams.patch_dist_limit_mult << "]\n";
            std::cerr << "\t [--angle-gens "; for(size_t vi=0;vi!=angle_gens.size();++vi)std::cerr<<angle_gens[vi]<<","; std::cerr << "]\n";
            std::cerr << "\t [--patch-pop-limit " << generatorParams.patch_population_limit << "]\n";
            std::cerr << "\t [--small-mode " << generatorParams.small_mode << "\t | 0: IGNORE, 1: RECEIVE_SIMILAR, 2: RECEIVE_ALL]\n";
            std::cerr << std::endl;

            if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
                return EXIT_FAILURE;
        }

        if ( boost::filesystem::is_directory(cloud_path) )
        {
            cloud_path += "/cloud.ply";
        }

        if ( !boost::filesystem::exists(cloud_path) )
        {
            std::cerr << "[" << __func__ << "]: " << "cloud file does not exist! " << cloud_path << std::endl;
            return EXIT_FAILURE;
        }
    } // ... parse input

    // Read desired angles
    if ( EXIT_SUCCESS == err )
    {
        processing::appendAnglesFromGenerators( generatorParams.angles, angle_gens, true );
    } //...read angles

    // Read points
    PointContainerT points;
    if ( EXIT_SUCCESS == err )
    {
        err = io::readPoints<PointPrimitiveT>( points, cloud_path );
        if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
    } //...read points

    std::vector<std::pair<int,int> > points_primitives;
    io::readAssociations( points_primitives, associations_path, NULL );
    for ( size_t i = 0; i != points.size(); ++i )
    {
        // store association in point
        points[i].setTag( PointPrimitiveT::GID, points_primitives[i].first );
    }

    // read primitives
    typedef std::map<int, typename PrimitiveContainerT::value_type> PrimitiveMapT;
    PrimitiveContainerT initial_primitives;
    PrimitiveMapT patches;
    {
        std::cout << "[" << __func__ << "]: " << "reading primitives from " << input_prims_path << "...";
        io::readPrimitives<PrimitiveT, typename PrimitiveContainerT::value_type>( initial_primitives, input_prims_path, &patches );
        std::cout << "reading primitives ok (#: " << initial_primitives.size() << ")\n";
    } //...read primitives

    //_____________________WORK_______________________
    //_______________________________________________

    // Generate
    //PrimitiveContainerT primitives;
    PrimitiveMapT primitives;
    if ( EXIT_SUCCESS == err )
    {
        err = CandidateGenerator::generate< MyPrimitivePrimitiveAngleFunctor, MyPointPrimitiveDistanceFunctor, PrimitiveT >
                                          ( primitives, patches/*initial_primitives*/, points, generatorParams.scale, generatorParams.angles, generatorParams );

        if ( err != EXIT_SUCCESS ) std::cerr << "[" << __func__ << "]: " << "generate exited with error! Code: " << err << std::endl;
    } //...generate

#if 0 // these shouldn't change here
    // Save point GID tags
    if ( EXIT_SUCCESS == err )
    {
        std::string assoc_path = boost::filesystem::path( cloud_path ).parent_path().string() + "/" + "points_primitives.csv";

        util::saveBackup( assoc_path );
        err = io::writeAssociations<PointPrimitiveT>( points, assoc_path );

        if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "saveBackup or writeAssociations exited with error! Code: " << err << std::endl;
        else                        std::cout << "[" << __func__ << "]: " << "wrote to " << assoc_path << std::endl;

    } //...save Associations
#endif

    // save primitives
    std::string o_path = boost::filesystem::path( cloud_path ).parent_path().string() + "/";
    if ( EXIT_SUCCESS == err )
    {
        std::string output_prims_path( o_path + "candidates.csv" );
        {
            int iteration = 0;
            iteration = util::parseIteration( input_prims_path ) + 1;
            std::stringstream ss;
            ss << o_path << "candidates_it" << iteration << ".csv";
            output_prims_path = ss.str();
        }

        util::saveBackup( output_prims_path );
        err = io::savePrimitives<PrimitiveT,typename PrimitiveContainerT::value_type::const_iterator>( /* what: */ primitives, /* where_to: */ output_prims_path );

        if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "saveBackup or savePrimitive exited with error! Code: " << err << std::endl;
        else                        std::cout << "[" << __func__ << "]: " << "wrote to " << output_prims_path << std::endl;
    } //...save primitives

    return err;
} // ...Solver::generate()

//! \brief      Step 3. Reads a formulated problem from path and runs qcqpcpp::OptProblem::optimize() on it.
//! \param argc Number of command line arguments.
//! \param argv Vector of command line arguments.
//! \return     Exit code. 0 == EXIT_SUCCESS.
int
Solver::solve( int    argc
             , char** argv )
{
    int                                   err           = EXIT_SUCCESS;

    bool                                  verbose       = false;
    enum SOLVER { MOSEK, BONMIN, GUROBI } solver        = MOSEK;
    std::string                           project_path  = "problem", solver_str = "bonmin";
    Scalar                                max_time      = 360;
    int                                   bmode         = 0; // Bonmin solver mode, B_Bb by default
    std::string                           rel_out_path  = ".";
    std::string                           x0_path       = "";

    // parse
    {
        bool valid_input = true;

        // verbose parsing
        if ( pcl::console::find_switch(argc,argv,"-v") || pcl::console::find_switch(argc,argv,"--verbose") )
            verbose = true;

        // solver parsing
        //std::string solver_str = "bonmin";
        valid_input &= pcl::console::parse_argument( argc, argv, "--solver", solver_str) >= 0;
        {
                 if ( !solver_str.compare("mosek")  ) solver = MOSEK;
            else if ( !solver_str.compare("bonmin") ) solver = BONMIN;
            else if ( !solver_str.compare("gurobi") ) solver = GUROBI;
            else
            {
                std::cerr << "[" << __func__ << "]: " << "Cannot parse solver " << solver_str << std::endl;
                valid_input = false;
            }
#           ifndef GF2_WITH_BONMIN
                 if ( solver == BONMIN )
                     throw new std::runtime_error("You have to specify a different solver, the project was not compiled with Bonmin enabled!");
#           endif // WITH_BONMIN
#           ifndef GF2_WITH_MOSEK
                 if ( solver == MOSEK )
                     throw new std::runtime_error("You have to specify a different solver, the project was not compiled with Mosek enabled!");
#           endif // WITH_MOSEK
#           ifndef GF2_WITH_GUROBI
                 if ( solver == GUROBI )
                     throw new std::runtime_error("You have to specify a different solver, the project was not compiled with GUROBI enabled!");
#           endif // WITH_GUROBI
        }
        // problem parsing
        pcl::console::parse_argument( argc, argv, "--problem", project_path );
        // parse time
        pcl::console::parse_argument( argc, argv, "--time", max_time );
        // parse bonmin solver mode
        pcl::console::parse_argument( argc, argv, "--bmode", bmode );
        pcl::console::parse_argument( argc, argv, "--rod"  , rel_out_path );

        // X0
        if (pcl::console::parse_argument( argc, argv, "--x0", x0_path ) >= 0)
        {
            if ( !boost::filesystem::exists(x0_path) )
            {
                std::cerr << "[" << __func__ << "]: " << "X0 path does not exist! " << x0_path << std::endl;
                valid_input = false;
            }
        }

        // usage print
        std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt\n"
                  << "\t--solver *" << solver_str << "* (mosek | bonmin | gurobi)\n"
                  << "\t--problem " << project_path << "\n"
                  << "\t[--time] " << max_time << "\n"
                  << "\t[--bmode *" << bmode << "*\n"
                         << "\t\t0 = B_BB, Bonmin\'s Branch-and-bound \n"
                         << "\t\t1 = B_OA, Bonmin\'s Outer Approximation Decomposition\n"
                         << "\t\t2 = B_QG, Bonmin's Quesada & Grossmann branch-and-cut\n"
                         << "\t\t3 = B_Hyb Bonmin's hybrid outer approximation\n"
                         << "\t\t4 = B_Ecp Bonmin's implemantation of ecp cuts based branch-and-cut a la FilMINT\n"
                         << "\t\t5 = B_IFP Bonmin's implemantation of iterated feasibility pump for MINLP]\n"
                  << "\t[--verbose] " << "\n"
                  << "\t[--rod " << rel_out_path << "]\t\t Relative output directory\n"
                  << "\t[--x0 " << x0_path << "]\t Path to starting point sparse matrix\n"
                  << "\t[--help, -h] "
                  << std::endl;

        // valid_input
        if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
        {
            std::cerr << "[" << __func__ << "]: " << "--solver is compulsory" << std::endl;
            err = EXIT_FAILURE;
        } //...if valid_input
    } //...parse

    // select solver
    typedef double OptScalar; // Mosek, and Bonmin uses double internally, so that's what we have to do...
    typedef qcqpcpp::OptProblem<OptScalar> OptProblemT;
    OptProblemT *p_problem = NULL;
    if ( EXIT_SUCCESS == err )
    {
        switch ( solver )
        {
#ifdef GF2_WITH_MOSEK
            case MOSEK:
                p_problem = new qcqpcpp::MosekOpt<OptScalar>( /* env: */ NULL );
                break;
#endif
            case BONMIN:
                p_problem = new qcqpcpp::BonminOpt<OptScalar>();
                break;

            default:
                std::cerr << "[" << __func__ << "]: " << "Unrecognized solver type, exiting" << std::endl;
                err = EXIT_FAILURE;
                break;
        } //...switch
    } //...select solver

    // problem.read()
    if ( EXIT_SUCCESS == err )
    {
        err += p_problem->read( project_path );
        if ( EXIT_SUCCESS != err )
            std::cerr << "[" << __func__ << "]: " << "Could not read problem, exiting" << std::endl;
    } //...problem.read()

    // problem.parametrize()
    {
        if ( max_time > 0 )
            p_problem->setTimeLimit( max_time );
        if ( solver == BONMIN )
        {
#           ifdef GF2_WITH_BONMIN
            static_cast<qcqpcpp::BonminOpt<OptScalar>*>(p_problem)->setAlgorithm( Bonmin::Algorithm(bmode) );
            OptProblemT::SparseMatrix x0;
            if ( !x0_path.empty() )
            {
                x0 = qcqpcpp::io::readSparseMatrix<OptScalar>( x0_path, 0 );
                static_cast<qcqpcpp::BonminOpt<OptScalar>*>(p_problem)->setStartingPoint( x0 );
            }

#           endif // WITH_BONMIN
        }
    }

    // problem.update()
    OptProblemT::ReturnType r = 0;
    if ( EXIT_SUCCESS == err )
    {
        // log
        if ( verbose ) { std::cout << "[" << __func__ << "]: " << "calling problem update..."; fflush(stdout); }

        // update
        r = p_problem->update();

#ifdef GF2_WITH_MOSEK
        // check output
        if ( r != MSK_RES_OK )
        {
            std::cerr << "[" << __func__ << "]: " << "ooo...update didn't work with code " << r << std::endl;
            err = EXIT_FAILURE;
        }
#endif

        // log
        if ( verbose ) { std::cout << "[" << __func__ << "]: " << "problem update finished\n"; fflush(stdout); }
    } //...problem.update()

    // problem.optimize()
    if ( EXIT_SUCCESS == err )
    {
        // optimize
        std::vector<OptScalar> x_out;
        if ( r == p_problem->getOkCode() )
        {
            // log
            if ( verbose ) { std::cout << "[" << __func__ << "]: " << "calling problem optimize...\n"; fflush(stdout); }

            // work
            r = p_problem->optimize( &x_out, OptProblemT::OBJ_SENSE::MINIMIZE );

            // check output
            if ( r != p_problem->getOkCode() )
            {
                std::cerr << "[" << __func__ << "]: " << "ooo...optimize didn't work with code " << r << std::endl; fflush(stderr);
                err = r;
            }
        } //...optimize

        // copy output
        std::vector<Scalar> scalar_x_out( x_out.size() );
        if ( EXIT_SUCCESS == err )
        {
            // copy
            std::copy( x_out.begin(), x_out.end(), scalar_x_out.begin() );

            // print energy
//            checkSolution( scalar_x_out
//                           , p_problem->getLinObjectivesMatrix().cast<Scalar>()
//                           , p_problem->getQuadraticObjectivesMatrix().cast<Scalar>()
//                           , p_problem->getLinConstraintsMatrix().cast<Scalar>()
//                           , weights );
        } //...copy output

        // dump
        {
            {
                std::string x_path = project_path + "/x.csv";
                OptProblemT::SparseMatrix sp_x( x_out.size(), 1 ); // output colvector
                for ( size_t i = 0; i != x_out.size(); ++i )
                {
                    if ( int(round(x_out[i])) > 0 )
                    {
                        sp_x.insert(i,0) = x_out[i];
                    }
                }
                qcqpcpp::io::writeSparseMatrix<OptScalar>( sp_x, x_path, 0 );
                std::cout << "[" << __func__ << "]: " << "wrote output to " << x_path << std::endl;
            }

            std::string candidates_path;
            if ( pcl::console::parse_argument( argc, argv, "--candidates", candidates_path ) >= 0 )
            {
                // read primitives
                PrimitiveContainerT prims;
                {
                    if ( verbose ) std::cout << "[" << __func__ << "]: " << "reading primitives from " << candidates_path << "...";
                    io::readPrimitives<PrimitiveT, typename PrimitiveContainerT::value_type>( prims, candidates_path );
                    if ( verbose ) std::cout << "reading primitives ok\n";
                } //...read primitives

                // save selected primitives
                PrimitiveContainerT out_prims( 1 );
                int prim_id = 0;
                for ( size_t l = 0; l != prims.size(); ++l )
                    for ( size_t l1 = 0; l1 != prims[l].size(); ++l1, ++prim_id )
                    {
                        if ( int(round(x_out[prim_id])) > 0 )
                        {
                            std::cout << "saving " << prims[l][l1].getTag(PrimitiveT::GID) << ", " << prims[l][l1].getTag(PrimitiveT::DIR_GID) << ", X: " << x_out[prim_id] << "\t, ";
                            prims[l][l1].setTag( PrimitiveT::CHOSEN, 1 );
                            out_prims.back().push_back( prims[l][l1] );
                        }
                    } // ... for l1
                std::cout << std::endl;

                std::string parent_path = boost::filesystem::path(candidates_path).parent_path().string();
                if ( parent_path.empty() )  parent_path = "./";
                else                        parent_path += "/";

                std::string out_prim_path = parent_path + rel_out_path + "/primitives." + solver_str + ".csv";
                {
                    int iteration = 0;
                    iteration = std::max(0,util::parseIteration(candidates_path) );
                    std::stringstream ss;
                    ss << parent_path + rel_out_path << "/primitives_it" << iteration << "." << solver_str << ".csv";
                    out_prim_path = ss.str();
                }

                util::saveBackup    ( out_prim_path );
                io::savePrimitives<PrimitiveT, PrimitiveContainerT::value_type::const_iterator>( out_prims, out_prim_path, /* verbose: */ true );
            } // if --candidates
            else
            {
                std::cout << "[" << __func__ << "]: " << "You didn't provide candidates, could not save primitives" << std::endl;
            } // it no --candidates
        }

    } //...problem.optimize()

    if ( p_problem ) { delete p_problem; p_problem = NULL; }

    return err;
}

//! \brief Unfinished function. Supposed to do GlobFit.
int
Solver::datafit( int    argc
               , char** argv )
{
    int                     err             = EXIT_SUCCESS;
    Scalar                  scale           = 0.05f;
    std::string             cloud_path      = "cloud.ply",
                            primitives_path = "candidates.csv",
                            associations_path = "";
    std::vector<Scalar>     angle_gens      = { Scalar(90.) };
    bool                    verbose         = false;

    // parse params
    {
        bool valid_input = true;
        valid_input &= pcl::console::parse_argument( argc, argv, "--scale"     , scale          ) >= 0;
        valid_input &= pcl::console::parse_argument( argc, argv, "--cloud"     , cloud_path     ) >= 0;
        valid_input &= pcl::console::parse_argument( argc, argv, "--primitives", primitives_path) >= 0;
        pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens );
        if ( pcl::console::find_switch(argc,argv,"-v") || pcl::console::find_switch(argc,argv,"--verbose") )
            verbose = true;

        if (    (pcl::console::parse_argument( argc, argv, "-a", associations_path) < 0)
             && (pcl::console::parse_argument( argc, argv, "--assoc", associations_path) < 0)
             && (!boost::filesystem::exists(associations_path)) )
        {
            std::cerr << "[" << __func__ << "]: " << "-a or --assoc is compulsory" << std::endl;
            valid_input = false;
        }

        std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt --gfit\n"
                  << "\t--scale " << scale << "\n"
                  << "\t--cloud " << cloud_path << "\n"
                  << "\t--primitives " << primitives_path << "\n"
                  << "\t--a,--assoc " << associations_path << "\n"
                  << "\t[--angle-gens "; for(size_t i=0;i!=angle_gens.size();++i)std::cerr<<angle_gens[i]<<",";std::cerr<< "]\n";
        std::cerr << std::endl;

        if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale, --cloud and --primitives are compulsory" << std::endl;
            return EXIT_FAILURE;
        }
    } // ... parse params
    // read points
    PointContainerT     points;
    {
        if ( verbose ) std::cout << "[" << __func__ << "]: " << "reading cloud from " << cloud_path << "...";
        io::readPoints<PointPrimitiveT>( points, cloud_path );
        std::vector<std::pair<int,int> > points_primitives;
        io::readAssociations( points_primitives, associations_path, NULL );
        for ( size_t i = 0; i != points.size(); ++i )
        {
            // store association in point
            points[i].setTag( PointPrimitiveT::GID, points_primitives[i].first );
        }
        if ( verbose ) std::cout << "reading cloud ok\n";
    } //...read points

    // read primitives
    typedef std::map<int, typename PrimitiveContainerT::value_type> PrimitiveMapT;
    PrimitiveContainerT prims;
    PrimitiveMapT       patches;
    {
        if ( verbose ) std::cout << "[" << __func__ << "]: " << "reading primitives from " << primitives_path << "...";
        io::readPrimitives<PrimitiveT, typename PrimitiveContainerT::value_type>( prims, primitives_path, &patches );
        if ( verbose ) std::cout << "reading primitives ok\n";
    } //...read primitives

    // Read desired angles
    std::vector<Scalar> angles;
    {
        processing::appendAnglesFromGenerators( angles, angle_gens, true );
    } // ... read angles

#if 0
    // WORK
    {
        typedef double                         OptScalar;
        typedef qcqpcpp::BonminOpt<OptScalar> OptProblemT;
        OptProblemT problem;

        const int Dims = 2; // 2: nx,ny, 3: +nz
        typedef typename PrimitiveContainerT::value_type::value_type PrimitiveT;        // LinePrimitive or PlanePrimitive
        typedef typename PointContainerT::value_type                 PointPrimitiveT;   // PointPrimitive to store oriented points
        typedef          std::pair<int,int>                          IntPair;
        typedef          std::map<IntPair, std::vector<int> >        PrimsVarsT;

        std::vector<std::pair<std::string,Scalar> >     starting_values;    // [var_id] = < var_name, x0 >
        PrimsVarsT                                      prims_vars;         // < <lid,lid1>, var_id >, associates primitives to variables

        char                                            name[255];          // variable name
        const Eigen::Matrix<Scalar,3,1>                 origin( Eigen::Matrix<Scalar,3,1>::Zero() ); // to calculate d

        typedef typename PrimitiveMapT::mapped_type InnerType;

        // add variables
        {
            int lid = 0;
            for ( PrimitiveMapT::const_iterator gid_it = patches.begin(); gid_it != patches.end(); ++gid_it, ++lid )
            {
                int lid1 = 0;
                for ( InnerType::const_iterator dir_it  = containers::valueOf<PrimitiveT>(gid_it).begin();
                                                dir_it != containers::valueOf<PrimitiveT>(gid_it).end();
                                              ++dir_it, ++lid1 )
                {
                    PrimitiveT const& prim = *dir_it;

                    const int gid     = prim.getTag( PrimitiveT::GID     ); //lid;
                    const int dir_gid = prim.getTag( PrimitiveT::DIR_GID ); //lid1;

                    Eigen::Matrix<Scalar,3,1> normal = prim.template normal<Scalar>();
                    std::cout << "line_" << gid << "_" << dir_gid << ".n = " << normal.transpose() << std::endl;

                    sprintf( name, "nx_%d_%d", gid, dir_gid );
                    prims_vars[ IntPair(lid,lid1) ].push_back( problem.addVariable(OptProblemT::BOUND::RANGE, -problem.getINF(), problem.getINF(), OptProblemT::VAR_TYPE::CONTINUOUS, OptProblemT::LINEAR, name) );
                    //vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name) );
                    starting_values.push_back( std::pair<std::string,Scalar>(name,normal(0)) );

                    sprintf( name, "ny_%d_%d", gid, dir_gid );
                    prims_vars[ IntPair(lid,lid1) ].push_back( problem.addVariable(OptProblemT::BOUND::RANGE, -problem.getINF(), problem.getINF(), OptProblemT::VAR_TYPE::CONTINUOUS, OptProblemT::LINEAR, name) );
                    //vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name) );
                    starting_values.push_back( std::pair<std::string,Scalar>(name,normal(1)) );

                    if ( Dims > 2 )
                    {
                        sprintf( name, "nz_%d_%d", gid, dir_gid );
                        //vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name) );
                        prims_vars[ IntPair(lid,lid1) ].push_back( problem.addVariable(OptProblemT::BOUND::RANGE, -problem.getINF(), problem.getINF(), OptProblemT::VAR_TYPE::CONTINUOUS, OptProblemT::LINEAR, name) );
                        starting_values.push_back( std::pair<std::string,Scalar>(name,normal(2)) );
                    }

                    sprintf( name, "d_%d_%d", gid, dir_gid );
                    Scalar d = Scalar(-1.) * prim.getDistance( origin );
                    std::cout << "line_" << gid << "_" << dir_gid << ".d = " << d << std::endl;
                    // vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name) );
                    prims_vars[ IntPair(lid,lid1) ].push_back( problem.addVariable(OptProblemT::BOUND::RANGE, -problem.getINF(), problem.getINF(), OptProblemT::VAR_TYPE::CONTINUOUS, OptProblemT::LINEAR, name) );
                    starting_values.push_back( std::pair<std::string,Scalar>(name,d) );
                }
            }
        }

        // add constraints: |normal|^2 = 1
        {
            int lid = 0;
            for ( PrimitiveMapT::const_iterator gid_it = patches.begin(); gid_it != patches.end(); ++gid_it, ++lid )
            {
                int lid1 = 0;
                for ( InnerType::const_iterator dir_it  = containers::valueOf<PrimitiveT>(gid_it).begin();
                                                dir_it != containers::valueOf<PrimitiveT>(gid_it).end();
                                              ++dir_it, ++lid1 )
                {
                    // sparse quadratic matrix
                    OptProblemT::SparseMatrix norm_constraint( problem.getVarCount(), problem.getVarCount() );

                    // 1 * nx * nx + 1 * ny * ny + 1 *  nz * nz
                    for ( int dim = 0; dim != Dims; ++dim )
                        norm_constraint.insert( prims_vars[IntPair(lid,lid1)][dim]
                                              , prims_vars[IntPair(lid,lid1)][dim] ) = Scalar( 1. );

                    // add constraint instance
                    problem.addConstraint  ( OptProblemT::BOUND::EQUAL, /* >= 1 */ Scalar(1.), /* <= 1 */ Scalar(1.), /* linear constraint coeffs: */ NULL );
                    // add quadratic coefficients
                    problem.addQConstraints( norm_constraint );
                }
            }
        }

        /// cost -> objective: minimize \sum_n \sum_p ((n_j . p_i) + d)^2  where point p_i is assigned to line with normal n_j
        if ( Dims == 3 ) throw new std::runtime_error("datafit unimplemented for 3D");

        {
            // assignments
            GidPidVectorMap populations;
            processing::getPopulations( populations, points );

            // for each line
            Scalar coeff = Scalar(0);
            int    lid   = 0;
            for ( PrimitiveMapT::const_iterator gid_it = patches.begin(); gid_it != patches.end(); ++gid_it, ++lid )
            {
                int lid1 = 0;
                int gid  = -2; // -1 is unset, -2 is unread
                for ( InnerType::const_iterator dir_it  = containers::valueOf<PrimitiveT>(gid_it).begin();
                                                dir_it != containers::valueOf<PrimitiveT>(gid_it).end();
                                              ++dir_it, ++lid1 )
                {
                    if ( gid < -1 )
                        gid = dir_it->getTag( PrimitiveT::GID );

                    // for each assigned point
                    for ( size_t pid_id = 0; pid_id != populations[gid].size(); ++pid_id )
                    {
                        const int pid = populations[gid][pid_id];

                        // skip, if ambiguous assignment
                        //if ( points_primitives[pid].size() != 1 )
                        //    continue;

                        // debug
                        if ( verbose ) std::cout << "[" << __func__ << "]: " << "adding pid " << pid << " -> " << "lines[" << lid << "][" << lid1 << "]" << std::endl;

                        // for each dimension: x,y,z
                        for ( int dim = 0; dim != Dims; ++dim )
                        {
                            // (p_x)^2 . (n_x)^2
                            coeff = points[pid].pos()( dim );                                                                                             // (p_x)
                            coeff *= coeff;                                                                                                               // (p_x)^2
                            problem.addQObjective( prims_vars[IntPair(lid,lid1)][dim], prims_vars[IntPair(lid,lid1)][dim], coeff );                       // (p_x)^2 . (n_x)^2

                            // debug
                            if ( verbose ) std::cout << "added qterm(pid: " << pid << "): "
                                                        << coeff << " * "
                                                        << problem.getVarName( prims_vars[IntPair(lid,lid1)][dim] ) << " * "
                                                        << problem.getVarName( prims_vars[IntPair(lid,lid1)][dim] )
                                                        << std::endl;

                            // 2 . p_x . n_x . d
                            coeff = Scalar(2.) * points[pid].pos()( dim );                                                                                // 2 . p_x
                            problem.addQObjective( /* n_x: */ prims_vars[IntPair(lid,lid1)][dim], /* d: */ prims_vars[IntPair(lid,lid1)][Dims], coeff );  // 2 . p_x . n_x . d

                            // debbug
                            if ( verbose ) std::cout << "added qterm(pid: " << pid << "): "
                                                        << coeff << " * "
                                                        << problem.getVarName( prims_vars[IntPair(lid,lid1)][dim] ) << " * "
                                                        << problem.getVarName( prims_vars[IntPair(lid,lid1)][Dims] )
                                                        << std::endl;
                        }

                        // d^2
                        coeff = Scalar(1);
                        problem.addQObjective( /* d: */ prims_vars[IntPair(lid,lid1)][Dims], /* d: */ prims_vars[IntPair(lid,lid1)][Dims], coeff );

                        // debug
                        if ( verbose ) std::cout << "added qterm(pid: " << pid << "): "
                                                    << coeff << " * "
                                                    << problem.getVarName( prims_vars[IntPair(lid,lid1)][Dims] ) << " * "
                                                    << problem.getVarName( prims_vars[IntPair(lid,lid1)][Dims] )
                                                    << std::endl;

#                       warning "[solver.h][datafit]TODO: derive this for pz,nz (Dims==3)"
                        // 2 . px . py . nx . ny
                        coeff = Scalar(2.) * points[pid].pos()(0) * points[pid].pos()(1);
                        problem.addQObjective( prims_vars[IntPair(lid,lid1)][0], prims_vars[IntPair(lid,lid1)][1], coeff );

                        // debug
                        if ( verbose ) std::cout << "added qterm(pid: " << pid << "): "
                                                    << coeff << " * "
                                                    << problem.getVarName( prims_vars[IntPair(lid,lid1)][0] ) << " * "
                                                    << problem.getVarName( prims_vars[IntPair(lid,lid1)][1] )
                                                    << std::endl;
                    } // for points
                } //...for directions
            } //...for patches
        } //...add objective

        // starting point
        {
            if ( starting_values.size() != problem.getVarCount() )
            {
                std::cerr << "[" << __func__ << "]: " << "starting_values.size() " << starting_values.size() << " != " << problem.getVarCount() << " problem.getVarCount()" << std::endl;
                return EXIT_FAILURE;
            }

            OptProblemT::SparseMatrix x0( problem.getVarCount(), 1 );
            for ( size_t i = 0; i != starting_values.size(); ++i )
            {
                x0.insert( i, 0 ) = starting_values[i].second;
            }
            problem.setStartingPoint( x0 );
        } //...starting values

        // save
        {
            problem.write( "./datafit_problem" );
        }

        // solve
        {
            OptProblemT::ReturnType r = problem.getOkCode();

            if ( problem.getOkCode() == r )
            {
                if ( verbose ) { std::cout << "[" << __func__ << "]: " << "calling problem optimize...\n"; fflush(stdout); }
                r = problem.update();
                if ( verbose ) { std::cout << "[" << __func__ << "]: " << "finished problem optimize...\n"; fflush(stdout); }
            }

            // optimize
            std::vector<OptScalar> x_out;
            if ( r == problem.getOkCode() )
            {
                // log
                if ( verbose ) { std::cout << "[" << __func__ << "]: " << "calling problem optimize...\n"; fflush(stdout); }

                // work
                r = problem.optimize( &x_out, OptProblemT::OBJ_SENSE::MINIMIZE );

                // check output
                if ( r != problem.getOkCode() )
                {
                    std::cerr << "[" << __func__ << "]: " << "ooo...optimize didn't work with code " << r << std::endl; fflush(stderr);
                    err = r;
                }
            } //...optimize
        }

    } //...work
#endif
    return err;
} // ...Solver::datafit()

//! \brief              Prints energy of solution in \p x using \p weights.
//! \param[in] x        A solution to calculate the energy of.
//! \param[in] weights  Problem weights used earlier. \todo Dump to disk together with solution.
Eigen::Matrix<Solver::Scalar,3,1>
Solver::checkSolution( std::vector<Scalar> const& x
                     , Solver::SparseMatrix const& linObj
                     , Solver::SparseMatrix const& Qo
                     , Solver::SparseMatrix const& /* A */
                     , Eigen::Matrix<Scalar,3,1> const& weights )
{
    Eigen::Matrix<Scalar,3,1> energy; energy.setZero();

    SparseMatrix complexity( x.size(), 1 );
    for ( size_t row = 0; row != x.size(); ++row )
        complexity.insert( row, 0 ) = weights(2);

    // X
    SparseMatrix mx( x.size(), 1 );
    for ( size_t i = 0; i != x.size(); ++i )
        mx.insert( i, 0 ) = x[i];

    SparseMatrix data   = linObj - complexity;

    // qo
    SparseMatrix e02 = mx.transpose() * linObj;
    std::cout << "[" << __func__ << "]: " << "qo * x = " << e02.coeffRef(0,0) << std::endl; fflush(stdout);

    // datacost
    SparseMatrix e0 = (mx.transpose() * data);
    energy(0) = e0.coeffRef(0,0);
    //std::cout << "[" << __func__ << "]: " << "data: " << energy(0) << std::endl; fflush(stdout);

    // Qo
    SparseMatrix e1 = mx.transpose() * Qo * mx;
    energy(1) = e1.coeffRef(0,0);
    //std::cout << "[" << __func__ << "]: " << "x' * Qo * x = pw = " << energy(1) << std::endl; fflush(stdout);

    // complexity
    SparseMatrix e2 = mx.transpose() * complexity;
    energy(2) = e2.coeffRef(0,0);
    //std::cout << "[" << __func__ << "]: " << "complx = " << energy(2) << std::endl; fflush(stdout);
    std::cout << "[" << __func__ << "]: " << std::setprecision(9) << energy(0) << " + " << energy(1) << " + " << energy(2) << " = " << energy.sum();
    std::cout                             << std::setprecision(9) << weights(0) << " * " << energy(0)/weights(0)
                                                                  << " + " << weights(1) << " * " << energy(1) / weights(1)
                                                                  << " + " << weights(2) << " * " << energy(2) / weights(2)
                                                                  << std::endl;

    return energy;
}

} // ... ns GF2

#endif // __GF2_SOLVER_H__
