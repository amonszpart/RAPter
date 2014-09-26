#ifndef GF2_SOLVER_HPP
#define GF2_SOLVER_HPP

#include "Eigen/Sparse"

#ifdef GF2_USE_PCL
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
//#   include "optimization/qp/gurobiOpt.h"
#endif

#include "globfit2/util/diskUtil.hpp"                 // saveBAckup
#include "globfit2/util/util.hpp"                     // timestamp2Str

#include "globfit2/io/io.h"
//#include "globfit2/optimization/candidateGenerator.h" // generate()
//#include "globfit2/optimization/energyFunctors.h"     // PointLineDistanceFunctor,
#include "globfit2/optimization/problemSetup.h"       // everyPatchNeedsDirection()

namespace GF2
{

//! \brief      Step 3. Reads a formulated problem from path and runs qcqpcpp::OptProblem::optimize() on it.
//! \param argc Number of command line arguments.
//! \param argv Vector of command line arguments.
//! \return     Exit code. 0 == EXIT_SUCCESS.
template <class _PrimitiveContainerT
         , class _InnerPrimitiveContainerT
         , class _PrimitiveT
         >
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
        valid_input &= (pcl::console::parse_argument( argc, argv, "--solver", solver_str) >= 0) || (pcl::console::parse_argument( argc, argv, "--solver3D", solver_str) >= 0);
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
    } //...problem.read()PrimitiveT

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
                _PrimitiveContainerT prims;
                {
                    if ( verbose ) std::cout << "[" << __func__ << "]: " << "reading primitives from " << candidates_path << "...";
                    io::readPrimitives<_PrimitiveT, _InnerPrimitiveContainerT>( prims, candidates_path );
                    if ( verbose ) std::cout << "reading primitives ok\n";
                } //...read primitives

                // save selected primitives
                _PrimitiveContainerT out_prims( 1 );
                int prim_id = 0;
                for ( size_t l = 0; l != prims.size(); ++l )
                    for ( size_t l1 = 0; l1 != prims[l].size(); ++l1 )
                    {
                        if ( prims[l][l1].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                        {
                            // copy small, keep for later iterations
                            out_prims.back().push_back( prims[l][l1] );
                        }
                        else
                        {
                            // copy to output, only, if chosen
                            if ( int(round(x_out[prim_id])) > 0 )
                            {
                                std::cout << "saving " << prims[l][l1].getTag(_PrimitiveT::GID) << ", " << prims[l][l1].getTag(_PrimitiveT::DIR_GID) << ", X: " << x_out[prim_id] << "\t, ";
                                prims[l][l1].setTag( _PrimitiveT::TAGS::STATUS, _PrimitiveT::STATUS_VALUES::ACTIVE );
                                out_prims.back().push_back( prims[l][l1] );
                            }

                            // increment non-small primitive ids
                            ++prim_id;
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
                io::savePrimitives<_PrimitiveT, typename _InnerPrimitiveContainerT::const_iterator>( out_prims, out_prim_path, /* verbose: */ true );
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
template < class _PrimitiveContainerT
         , class _InnerPrimitiveContainerT
         , class _PrimitiveT
>
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
                  << "\t--no-paral\n"
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
    typedef std::map<int, _InnerPrimitiveContainerT> PrimitiveMapT;
    _PrimitiveContainerT prims;
    PrimitiveMapT       patches;
    {
        if ( verbose ) std::cout << "[" << __func__ << "]: " << "reading primitives from " << primitives_path << "...";
        io::readPrimitives<_PrimitiveT, _InnerPrimitiveContainerT>( prims, primitives_path, &patches );
        if ( verbose ) std::cout << "reading primitives ok\n";
    } //...read primitives

    // Read desired angles
    bool no_paral = pcl::console::find_switch( argc, argv, "--no-paral");
    std::vector<Scalar> angles;
    {
        processing::appendAnglesFromGenerators( angles, angle_gens, no_paral, true );
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

//! \brief              Prints energy of solution in \p x using \p weights. Unused for now.
//! \param[in] x        A solution to calculate the energy of.
//! \param[in] weights  Problem weights used earlier. \todo Dump to disk together with solution.
Eigen::Matrix<GF2::Scalar,3,1>
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
} //...checkSolution

} // ... ns GF2

#endif // GF2_SOLVER_HPP
