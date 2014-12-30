#ifndef GF2_PROBLEMSETUP_HPP
#define GF2_PROBLEMSETUP_HPP

#include <iostream>                               // cout, cerr, endl
#include "Eigen/Dense"
#include <vector>
#include <map>
#include <set> // edgelist

#if GF2_USE_PCL
#   include "pcl/console/parse.h" // pcl::console::parse_argument
#endif

#include "qcqpcpp/optProblem.h"                   // OptProblem

#include "globfit2/optimization/energyFunctors.h" // AbstractPrimitivePrimitiveEnergyFunctor,
#include "globfit2/parameters.h"                  // ProblemSetupParams
#include "globfit2/processing/util.hpp"           // getPopulation()
#include "globfit2/io/io.h"                       // readPrimitives(), readPoints()

#include "globfit2/processing/graph.hpp"

namespace GF2 {

//____________________________________________class ProblemSetup __________________________________________________

template < class _PrimitiveContainerT
         , class _PointContainerT
         , class _PrimitiveT
         , class _PointPrimitiveT
         , class _FiniteFiniteDistFunctor
         > int
ProblemSetup::formulateCli( int    argc
                          , char** argv )
{
    bool verbose = false;
    typedef          MyPointPrimitiveDistanceFunctor               _PointPrimitiveDistanceFunctor;
    //typedef          MyPointPrimitiveDistanceFunctor               _PointPrimitiveDistanceFunctor;

    typedef typename _PointPrimitiveT::Scalar                      Scalar;

    GF2::ProblemSetupParams<Scalar> params;

    std::string               cloud_path         = "cloud.ply",
                              candidates_path    = "candidates.csv",
                              assoc_path         = "points_primitives.csv";
    std::vector<Scalar>       angle_gens         = { 90. };
    int                       srand_val          = 123456;
    std::string               cost_string        = "sqrt";
    std::string               problem_rel_path   = "problem";
    std::string               data_cost_mode_str = "assoc";
    std::string               constr_mode_str    = "hybrid";
    std::string               energy_path        = "energy.csv";
    bool                      calc_energy        = false; // instead of writing the problem, calculate the energy of selecting all input lines.
    // parse params
    {
        bool valid_input = true;
        verbose = pcl::console::find_switch( argc, argv, "--verbose") || pcl::console::find_switch( argc, argv, "-v" );

        valid_input &= pcl::console::parse_argument( argc, argv, "--scale"     , params.scale   ) >= 0;
        // cloud
        pcl::console::parse_argument( argc, argv, "--cloud"     , cloud_path     );
        valid_input &= boost::filesystem::exists( cloud_path );

        valid_input &= pcl::console::parse_argument( argc, argv, "--candidates", candidates_path) >= 0;
        pcl::console::parse_argument( argc, argv, "--unary", params.weights(0) );
        pcl::console::parse_argument( argc, argv, "--pw"   , params.weights(1) );
        pcl::console::parse_argument( argc, argv, "--cmp"  , params.weights(2) );
        pcl::console::parse_argument( argc, argv, "--cost-fn", cost_string );
        pcl::console::parse_argument( argc, argv, "--srand", srand_val );
        pcl::console::parse_argument( argc, argv, "--rod"  , problem_rel_path );
        pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens );
        pcl::console::parse_argument( argc, argv, "--dir-bias", params.dir_id_bias );
        if ( (pcl::console::parse_argument( argc, argv, "--assoc", assoc_path) < 0) && pcl::console::parse_argument( argc, argv, "-a", assoc_path ) < 0 )
        {
            std::cerr << "[" << __func__ << "]: " << "associations (points_primitives.csv) is compulsory!" << std::endl;
            valid_input = false;
        }

        // data_cost_mode
        {
            pcl::console::parse_argument( argc, argv, "--data-mode", data_cost_mode_str );
            /**/ if ( !data_cost_mode_str.compare( "assoc"   ) ) params.data_cost_mode = ProblemSetupParams<Scalar>::ASSOC_BASED; // default
            else if ( !data_cost_mode_str.compare( "band"    ) ) params.data_cost_mode = ProblemSetupParams<Scalar>::BAND_BASED;  // banded
            else if ( !data_cost_mode_str.compare( "instance") ) params.data_cost_mode = ProblemSetupParams<Scalar>::INSTANCE_BASED;
            else std::cerr << "[" << __func__ << "]: " << "could NOT parse data mode, assuming " << (int)params.data_cost_mode << std::endl;
        }

        // constr_mode
        {
            pcl::console::parse_argument( argc, argv, "--constr-mode", constr_mode_str );
            if ( verbose ) std::cout << "parsed constr_mode_str: " << constr_mode_str << std::endl;
            /**/ if ( !constr_mode_str.compare( "patch"  ) ) params.constr_mode = ProblemSetupParams<Scalar>::PATCH_WISE; // after first iteration
            else if ( !constr_mode_str.compare( "point"  ) ) params.constr_mode = ProblemSetupParams<Scalar>::POINT_WISE; // banded
            else if ( !constr_mode_str.compare( "hybrid" ) ) params.constr_mode = ProblemSetupParams<Scalar>::HYBRID; // hybrid, default
            else
            {
                std::cerr << "[" << __func__ << "]: " << "could NOT parse constraint mode (" << constr_mode_str << "), assuming " << (int)params.constr_mode << std::endl;
                throw new std::runtime_error("Could not parse constr-mode");
            }
        }

        // pop_limit
        {
            pcl::console::parse_argument( argc, argv, "--patch-pop-limit", params.patch_population_limit );
        }

        // energy
        {
            calc_energy = pcl::console::find_switch( argc, argv, "--energy" );
        }

        // freq_weight
        pcl::console::parse_argument( argc, argv, "--freq-weight", params.freq_weight );


        if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") || verbose )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale, --cloud and --candidates are compulsory" << std::endl;
            std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt --formulate\n"
                      << " --scale " << params.scale << "\n"
                      << " --cloud " << cloud_path << "\n"
                      << " --candidates " << candidates_path << "\n"
                      << " -a,--assoc " << assoc_path << "\n"
                      << " [--energy " << (calc_energy?"yes":"no") << "]\n"
                      << " [--angle-gens "; for(size_t vi=0;vi!=angle_gens.size();++vi)std::cerr<<angle_gens[vi]<<","; std::cerr << "]\n";
            std::cerr << " [--dir-bias " << params.dir_id_bias << "]\n"
                      << " [--unary " << params.weights(0) << "]\n"
                      << " [--pw " << params.weights(1) << "]\n"
                      << " [--cmp " << params.weights(2) << "]\n"
                      << " [--cost-fn *" << cost_string << "* (spatsqrt | cexp | sqrt)" << "]\n"
                      << " [--data-mode *" << (int)params.data_cost_mode << "* (assoc | band | instance) ]\n"
                      << " [--constr-mode *" << (int)params.constr_mode << "* (patch | point | hybrid ) ]\n"
                      << " [--srand " << srand_val << "]\n"
                      << " [--rod " << problem_rel_path << "]\t\tRelative output path of the output matrix files\n"
                      << " [--patch-pop-limit " << params.patch_population_limit << "]\n"
                      << " [--freq-weight" << params.freq_weight << "]\n"
                      << " [--energy-out" << energy_path << "]\n"
                      << " [--no-paral\n]"
                      << std::endl;
            if ( !verbose )
                return EXIT_FAILURE;
        }
    } // ... parse params

    // Read desired angles
    bool no_paral = pcl::console::find_switch( argc, argv, "--no-paral");
    {
        processing::appendAnglesFromGenerators( params.angles, angle_gens, no_paral, verbose );
    } //...read angles

    // read points
    _PointContainerT     points;
    {
        if ( verbose ) std::cout << "[" << __func__ << "]: " << "reading cloud from " << cloud_path << "...";
        io::readPoints<_PointPrimitiveT>( points, cloud_path );
        if ( verbose ) std::cout << "reading cloud ok\n";
    } //...read points

    // read primitives
    _PrimitiveContainerT prims;
    {
        if ( verbose ) std::cout << "[" << __func__ << "]: " << "reading primitives from " << candidates_path << "...";
        io::readPrimitives<_PrimitiveT, typename _PrimitiveContainerT::value_type>( prims, candidates_path );
        if ( verbose ) std::cout << "reading primitives ok\n";
    } //...read primitives

    // read point-primitive associations
    std::vector<std::pair<int,int> > points_primitives;
    io::readAssociations( points_primitives, assoc_path, NULL );
    for ( size_t i = 0; i != points.size(); ++i )
    {
        if ( i > points_primitives.size() )
        {
            std::cerr << "more points than associations..." << std::endl;
            return EXIT_FAILURE;
        }

        points[i].setTag( _PointPrimitiveT::GID, points_primitives[i].first );

        if ( points[i].getTag(_PointPrimitiveT::GID) == -1 )
            std::cerr << "[" << __func__ << "]: " << "point assigned to patch with id -1" << std::endl;
    } //...read associations

    //AbstractPrimitivePrimitiveEnergyFunctor<Scalar,_PrimitiveT> *primPrimDistFunctor = NULL;
    SpatialSqrtPrimitivePrimitiveEnergyFunctor<_FiniteFiniteDistFunctor, _PointContainerT, Scalar,_PrimitiveT> *primPrimDistFunctor = NULL;
    // parse cost function
    {
        if ( !cost_string.compare("spatsqrt") )
        {
            primPrimDistFunctor = new SpatialSqrtPrimitivePrimitiveEnergyFunctor<_FiniteFiniteDistFunctor, _PointContainerT, Scalar,_PrimitiveT>
                    ( params.angles, points, params.scale );
            primPrimDistFunctor->_verbose = verbose;
        }
        else
        {
            std::cerr << "[" << __func__ << "]: " << "unrecognized primitive-primitive cost function: " << cost_string << " you probably need --cost-fn \"spatsqrt\""<< std::endl;
            throw new std::runtime_error( "Unrecognized primitive-primitive cost function" );
        }
//             if ( !cost_string.compare("cexp") )    primPrimDistFunctor = new CExpPrimitivePrimitiveEnergyFunctor<Scalar,_PrimitiveT>( params.angles );
//        else if ( !cost_string.compare("sqrt") )    primPrimDistFunctor = new SqrtPrimitivePrimitiveEnergyFunctor<Scalar,_PrimitiveT>( params.angles );
//        else                                        std::cerr << "[" << __func__ << "]: " << "Could not parse cost functor input..." << std::endl;
    } //...parse cost function

    // WORK
    problemSetup::OptProblemT problem;
    AnglesT angle_gens_in_rad;
    for ( AnglesT::const_iterator angle_it = angle_gens.begin(); angle_it != angle_gens.end(); ++angle_it )
        angle_gens_in_rad.push_back( *angle_it * M_PI / 180. );
    int err = formulate<_PointPrimitiveDistanceFunctor>( problem
                                                       , prims
                                                       , points
                                                       , params.constr_mode
                                                       , params.data_cost_mode
                                                       , params.scale
                                                       , params.weights
                                                       , primPrimDistFunctor
                                                       , angle_gens_in_rad
                                                       , params.patch_population_limit
                                                       , params.dir_id_bias
                                                       , !calc_energy && verbose
                                                       , params.freq_weight );

    // dump. default output: ./problem/*.csv; change by --rod
    if ( EXIT_SUCCESS == err )
    {
        std::string parent_path = boost::filesystem::path(candidates_path).parent_path().string();
        if ( !parent_path.empty() )     parent_path += "/";
        else                            parent_path =  ".";

        if ( !calc_energy )
        {
            std::string problem_path = parent_path + "/" + problem_rel_path;
            problem.write( problem_path );
        }
        else
        {
            problemSetup::OptProblemT::VectorX x( problem.getVarCount(), 1 );
            x.setOnes();
            Scalar dataPlusCmplx = (x.transpose() * problem.getLinObjectivesMatrix()).coeff(0);
            Scalar complexityC   = x.sum() * params.weights(2);
            Scalar dataC         = dataPlusCmplx - complexityC;
            Scalar pairwiseC     = (x.transpose() * problem.getQuadraticObjectivesMatrix() * x );
            std::cout << "E = " << dataC + pairwiseC + complexityC << " = "
                      << dataC << " (data) + " << pairwiseC << "(pw) + " << complexityC << " (cmplx)"
                      << " = "
                      << params.weights(0) << " * " << dataC / params.weights(0)
                      << " + " << params.weights(1) << " * " << pairwiseC / params.weights(1)
                      << " + " << params.weights(2) << " * " << complexityC / params.weights(2)
                      << std::endl;
            std::ofstream fenergy( parent_path + "/" + energy_path, std::ofstream::out | std::ofstream::app );
            fenergy << dataC + pairwiseC + complexityC << "," << dataC << "," << pairwiseC << "," << complexityC << std::endl;
            fenergy.close();
        }
    } //...dump

    // cleanup
    if ( primPrimDistFunctor ) { delete primPrimDistFunctor; primPrimDistFunctor = NULL; }

    return err;
} //...ProblemSetup::formulateCli()

namespace problemSetup
{
    /*! \brief tuple<3> to temporarily store graph edges until we know how large the graph will be
     */
    template <typename _Scalar>
    struct EdgeT
    {
        int _v0, _v1;
        _Scalar _w;
        EdgeT( int v0, int v1, _Scalar w ) : _v0(v0), _v1(v1), _w(w) {}

        bool operator<( EdgeT<_Scalar> const& other ) const
        {
            if ( _v0 < other._v0 ) return true;
            else if ( _v1 < other._v1 ) return true;
            else return _w < other._w;
        }
    };

} //...ns problemSetup

template < class _PointPrimitiveDistanceFunctor
         , class _PrimitiveContainerT
         , class _PointContainerT
         , class _PrimPrimDistFunctorT
         , class _PrimitiveT
         , class _PointPrimitiveT
         , typename _Scalar
         > int
ProblemSetup::formulate( problemSetup::OptProblemT                                          & problem
                       , _PrimitiveContainerT                                          const& prims
                       , _PointContainerT                                              const& points
                       , typename ProblemSetupParams<_Scalar>::CONSTR_MODE             const  constr_mode
                       , typename ProblemSetupParams<_Scalar>::DATA_COST_MODE          const  data_cost_mode
                       , _Scalar                                                       const  scale
                       , Eigen::Matrix<_Scalar,-1,1>                                   const& weights
                       , _PrimPrimDistFunctorT                                       * const& primPrimDistFunctor // AbstractPrimitivePrimitiveEnergyFunctor<_Scalar,_PrimitiveT>
                       , AnglesT                                                       const& angle_gens_in_rad
                       , int                                                           const  patch_pop_limit
                       , _Scalar                                                       const  dir_id_bias
                       , int                                                           const  verbose
                       , _Scalar                                                       const  freq_weight /* = 0. */
        )
{
    using problemSetup::OptProblemT;

    typedef Graph< _Scalar, typename MyGraphConfig<_Scalar>::UndirectedGraph > GraphT;
    typedef problemSetup::EdgeT<_Scalar> EdgeT;

    // work - formulate problem
    int err = EXIT_SUCCESS;

    // log
    if ( verbose ) { std::cout << "[" << __func__ << "]: " << "formulating problem...\n"; fflush(stdout); }

    typedef std::pair<int,int> IntPair;
    /*   */ std::map <IntPair,int> lids_varids;
    std::set<int> chosen_varids;

    AnglesT angles = primPrimDistFunctor->getAngles();
    if ( angle_gens_in_rad.end() != std::find_if( angle_gens_in_rad.begin(), angle_gens_in_rad.end(), [](Scalar const& angle) { return angle > 2. * M_PI; } ) )
    {
        std::cerr << "[" << __func__ << "]: " << "angle_gens need to be in rad, are you sure?" << std::endl;
        std::cerr<<"angle_gens_in_rad:";for(size_t vi=0;vi!=angle_gens_in_rad.size();++vi)std::cerr<<angle_gens_in_rad[vi]<<" ";std::cerr << "\n";
        throw new std::runtime_error("angle_gens need to be in rad, are you sure");
    }


    // ____________________________________________________
    // Variables - Add variables to problem
    {
        char name[16];
        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                // ignore smalls
                if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                    continue;

                // add var
                int gid     = prims[lid][lid1].getTag( _PrimitiveT::TAGS::GID );
                int dir_gid = prims[lid][lid1].getTag( _PrimitiveT::TAGS::DIR_GID );
                sprintf( name, "x_%d_%d", gid, dir_gid );

                // store var_id for later, add binary variable
                const int var_id = problem.addVariable( OptProblemT::BOUND::RANGE, 0.0, 1.0, OptProblemT::VAR_TYPE::INTEGER
                                                        , OptProblemT::LINEARITY::LINEAR, name ); // changed to nonlinear by Aron on 29.12.2014
                lids_varids[ IntPair(lid,lid1) ] = var_id;

                // save for initial starting point
                if ( prims[lid][lid1].getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::ACTIVE )
                    chosen_varids.insert( var_id );
            }
        }
    } // ... variables

    // ____________________________________________________
    // Lin constraints
    if ( EXIT_SUCCESS == err )
    {
        switch ( constr_mode )
        {
            case ProblemSetupParams<_Scalar>::CONSTR_MODE::PATCH_WISE:
                err = problemSetup::everyPatchNeedsDirectionConstraint<_PointPrimitiveDistanceFunctor, _PrimitiveT, _PointPrimitiveT>
                        ( problem, prims, points, lids_varids, weights, scale, verbose );
                break;

            case ProblemSetupParams<_Scalar>::CONSTR_MODE::POINT_WISE:
                err = problemSetup::everyPointNeedsPatchConstraint<_PointPrimitiveDistanceFunctor, _PrimitiveT, _PointPrimitiveT>
                        ( problem, prims, points, lids_varids, weights, scale );
                break;
            case ProblemSetupParams<_Scalar>::CONSTR_MODE::HYBRID:
                err = problemSetup::largePatchesNeedDirectionConstraint<_PointPrimitiveDistanceFunctor, _PrimitiveT, _PointPrimitiveT>
                        ( problem, prims, points, lids_varids, weights, scale, patch_pop_limit, verbose );
                break;
        } //...switch constr_mode

        // Error check - Lin constraints
        if ( EXIT_SUCCESS != err )
        {
            std::cerr << "[" << __func__ << "]: " << "lin constraints setup returned error " << err << std::endl;
            return err;
        }
    } //...Lin constraints

    // ____________________________________________________
    // Unary cost -> lin objective
    if ( EXIT_SUCCESS == err )
    {
        switch ( data_cost_mode )
        {
            case ProblemSetupParams<_Scalar>::DATA_COST_MODE::ASSOC_BASED:
                err = problemSetup::associationBasedDataCost<_PointPrimitiveDistanceFunctor, _PrimitiveT, _PointPrimitiveT>
                        ( problem, prims, points, lids_varids, weights, scale, freq_weight, verbose );
                break;

            case ProblemSetupParams<_Scalar>::DATA_COST_MODE::INSTANCE_BASED:
                std::cerr << "::INSTANCE_BASED not implemented" << std::endl;
                throw new std::runtime_error("::INSTANCE_BASED not implemented");
//                err = problemSetup::instanceBasedDataCost<_PointPrimitiveDistanceFunctor, _PrimitiveT, _PointPrimitiveT>
//                        ( problem, prims, points, lids_varids, weights, scale );
                break;

            case ProblemSetupParams<_Scalar>::DATA_COST_MODE::BAND_BASED:
            default:
                err = problemSetup::bandBasedDataCost<_PointPrimitiveDistanceFunctor, _PrimitiveT, _PointPrimitiveT>
                        ( problem, prims, points, lids_varids, weights, scale );
                break;
        } //...switch data_cost_mode
    } //...unary cost

    // Error check - Unary costs
    if ( EXIT_SUCCESS != err )
    {
        std::cerr << "[" << __func__ << "]: " << "data cost setup returned error " << err << std::endl;
        return err;
    }

    // ____________________________________________________
    // Pairwise cost -> quad objective
    if ( EXIT_SUCCESS == err )
    {
        typedef std::vector< Eigen::Matrix<_Scalar,3,1> > ExtremaT;
        typedef std::pair<int,int> LidLid;
        typedef std::map< LidLid, ExtremaT > ExtremaMapT;
        ExtremaMapT extremas;
        typedef OptProblemT::Scalar ProblemScalar;

        GraphT::testGraph();
        std::set< EdgeT > edgesList;

        GidPidVectorMap populations;
        processing::getPopulations( populations, points );

        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                    continue;

                const int gid = prims[lid][lid1].getTag( _PrimitiveT::TAGS::GID );
                const int did = prims[lid][lid1].getTag( _PrimitiveT::TAGS::DIR_GID );

                // extremas key
                LidLid lidLid1( lid, lid1 );
                // find/calculate extrema
                typename ExtremaMapT::const_iterator it = extremas.find(lidLid1);
                if ( it == extremas.end() )
                {
                    prims[lid][lid1].template getExtent<_PointPrimitiveT>( extremas[lidLid1]
                                                       , points
                                                       , scale
                                                       , &(populations[gid]) );
                    it = extremas.find( lidLid1 );
                }

                for ( size_t lidOth = 0; lidOth != prims.size(); ++lidOth )
                {
                    for ( size_t lid1Oth = 0; lid1Oth != prims[lidOth].size(); ++lid1Oth )
                    {
                        if ( prims[lidOth][lid1Oth].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                            continue;

                        // skip same line, that's always zero
                        if ( (lid == lidOth) && (lid1 == lid1Oth) ) continue;

                        const int gidOth = prims[lidOth][lid1Oth].getTag( _PrimitiveT::TAGS::GID );
                        const int didOth = prims[lidOth][lid1Oth].getTag( _PrimitiveT::TAGS::DIR_GID );

                        // extremas key
                        LidLid lidLid1Oth( lidOth, lid1Oth );
                        // find/calculate extrema
                        typename ExtremaMapT::const_iterator oit = extremas.find( lidLid1Oth );
                        if ( oit == extremas.end() )
                        {
                            prims[lidOth][lid1Oth].template getExtent<_PointPrimitiveT>( extremas[lidLid1Oth]
                                                               , points
                                                               , scale
                                                               , &(populations[gidOth]) );
                            oit = extremas.find( lidLid1Oth );
                        }

                        // SqrtAngle:
                        _Scalar score = primPrimDistFunctor->eval( prims[lid][lid1], (*it).second
                                                                 , prims[lidOth][lid1Oth], (*oit).second
                                                                 , /* allowedAngles.size() ? allowedAngles[didOth] : */ angles );

                        // should be deprecated:
                        if ( prims[lid][lid1].getTag(_PrimitiveT::DIR_GID) != prims[lidOth][lid1Oth].getTag(_PrimitiveT::DIR_GID) )
                            score += dir_id_bias;

                        // multiply by pairwise weight
                        score *= weights(1);

                        // add
                        const int varId0 = lids_varids.at( IntPair(lid,lid1) );
                        const int varId1 = lids_varids.at( IntPair(lidOth,lid1Oth) );
                        problem.addQObjective( varId0, varId1, score );

                        // coupling
                        {
                            _Scalar invDist = primPrimDistFunctor->evalSpatial( prims[lid][lid1], (*it).second
                                                                              , prims[lidOth][lid1Oth], (*oit).second );
                            if ( invDist > _Scalar(0.) && (did == didOth) )
                            {
                                edgesList.insert( EdgeT(varId0, varId1, invDist) );
                            }
                        }
                    } // ... olid1
                } // ... olid
            } // ... lid1
        } // ... lid

        {
            typedef std::vector<int> ComponentListT;
            GraphT graph( lids_varids.size() );
            for ( auto it = edgesList.begin(); it != edgesList.end(); ++it )
            {
                graph.addEdge( it->_v0, it->_v1, /* not used right now: */ it->_w );
            }

            for ( int i = 0; i != lids_varids.size(); ++i )
            {
                if ( !problem.getVarName(i).empty() )
                    graph.addVertexName( i, problem.getVarName(i) );
            }

            {
                std::ofstream f;
                f.open( "components.gv" );
                f << "graph {\n";
                ComponentListT components;
                std::map<int,int> compSizes;
                graph.getComponents( components, &compSizes );

                std::map< int, std::vector<int> > clusters; // [ cluster0: [v0, v10,...], cluster1: [v3, v5, ...], ... ]
                for ( int varId = 0; varId != components.size(); ++varId )
                {
                    if ( compSizes[ components[varId] ] < 2 )
                        continue;

                    if ( problem.getVarName(varId).empty() )
                        f << varId;
                    else
                        f << problem.getVarName(varId);

                    f << " -- " << components[varId] << std::endl;

                    // store
                    clusters[ components[varId] ].push_back( varId );
                }

                f << "}" << std::endl;
                f.close();

                if ( clusters.size() )
                {
                    //system( "dot -Tpng -o comps.png components.gv && (eog comps.png &)" );

                    // work
                    typename OptProblemT::SparseMatrix cluster_constraint( 1, problem.getVarCount()+1 );
                    char name[255];
                    int clusterId = 0;
                    for ( auto it = clusters.begin(); it != clusters.end(); ++it, ++clusterId )
                    {
                        sprintf(name,"cl_%d", clusterId );
                        int varid = problem.addVariable( OptProblemT::BOUND::RANGE, 0.0, 1.0, OptProblemT::VAR_TYPE::INTEGER
                                                         , OptProblemT::LINEARITY::LINEAR, name );

                        for ( int i = 0; i != it->second.size(); ++i )
                        {
                            cluster_constraint.insert( 0, it->second[i] ) = 1;
                        }
                        cluster_constraint.insert( 0, varid ) = -(int)it->second.size();
                        // k * X_cluster_l <= A( l, : ) * X <= INF
                        problem.addConstraint( OptProblemT::BOUND::GREATER_EQ
                                               , /* lower_limit: */ 0 // k * X_cluster_l
                                               , /* upper_limit: */ problem.getINF()
                                               , &cluster_constraint );

                        chosen_varids.insert( varid );

                        // extra!
                        //problem.addLinObjective( varid, -1 );
                    }
                }
            }

            graph.draw( "graph.gv" );
            //system( "dot -Tpng -o graph.png graph.gv && (eog graph.png &)" );
        }

    } //...pairwise cost

    // ____________________________________________________
    // Initial solution
    {
        if ( chosen_varids.size() )
        {
            OptProblemT::SparseMatrix x0( problem.getVarCount(), 1 );
            for ( std::set<int>::const_iterator it = chosen_varids.begin(); it != chosen_varids.end(); ++it )
            {
                x0.insert( *it, 0 ) = 1;
            }
            problem.setStartingPoint( x0 );
        }
    } //...Initial solution

    // log
    if ( (EXIT_SUCCESS == err) && verbose )
    {
        problem.printProblem();
        std::cout << "[" << __func__ << "]: " << "formulating problem finished\n";
        fflush(stdout);
    } //...log

    return err;
} //...ProblemSetup::formulate()


//____________________________________________namespace problemSetup __________________________________________________

namespace problemSetup {

    //____________________________________________Constraints__________________________________________________

    //! \brief              Adds constraints to \p problem so, that each point has at least one line in \p scale radius that is selected.
    //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
    template < class _PointPrimitiveDistanceFunctor
             , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
             , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
             , typename _Scalar
             , class _OptProblemT
             , class _PrimitiveContainerT
             , class _PointContainerT
             , class _AssocT
             , class _WeightsT
             >
    static inline int
    everyPointNeedsPatchConstraint( _OptProblemT              & problem
                                  , _PrimitiveContainerT const& prims
                                  , _PointContainerT     const& points
                                  , _AssocT              const& lids_varids
                                  , _WeightsT            const& /*weights*/
                                  , _Scalar              const  scale )
    {
        typedef typename _AssocT::key_type      IntPair;
        typedef typename _OptProblemT::Scalar   ProblemScalar;

        int err = EXIT_SUCCESS;

        std::set< std::vector<ProblemScalar> > uniqueA;
        // each point needs at least one line
        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            // work
            std::vector<ProblemScalar> coeffs( problem.getVarCount(), 0 );
            bool non_zero = false;
            for ( size_t lid = 0; lid != prims.size(); ++lid )
                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                {
                    if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                        continue;

                    _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );
                    if ( dist < scale )
                    {
                        coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = 1.0;
                        non_zero = true;
                    }
                }

            if ( non_zero )
            {
                unsigned prev_size = uniqueA.size();
                uniqueA.insert( coeffs );
                if ( uniqueA.size() > prev_size )
                    problem.addLinConstraint( _OptProblemT::BOUND::GREATER_EQ, 1, problem.getINF(), coeffs );
            }
        } //... for each point

        return err;
    } //...everyPointNeedsPatchConstraint

    //! \brief              Adds constraints to \p problem so, that each patch (prims[i] that have the same _PrimitiveT::GID) has at least one member j (prims[i][j]) selected.
    //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
    template < class _PointPrimitiveDistanceFunctor
             , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
             , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
             , typename _Scalar
             , class _OptProblemT
             , class _PrimitiveContainerT
             , class _PointContainerT
             , class _AssocT
             , class _WeightsT
             >
    static inline int
    everyPatchNeedsDirectionConstraint( _OptProblemT              & problem
                                      , _PrimitiveContainerT const& prims
                                      , _PointContainerT     const& /*points*/
                                      , _AssocT              const& lids_varids
                                      , _WeightsT            const& /*weights*/
                                      , _Scalar              const  /*scale*/
                                      , bool                 const  verbose     /* = false */ )
    {
        typedef typename _AssocT::key_type IntPair;

        int err = EXIT_SUCCESS;

        // one direction / patch needs to be choosen
        std::vector< double              > coeffs ( problem.getVarCount(), 0 ); // constraint coefficient line in A
        std::set   < std::vector<double> > uniqueA;                             // to ensure unique constraints
        // for all patches
        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            bool non_zero_line = false;

            // init row to 0
            std::fill( coeffs.begin(), coeffs.end(), 0 );

            // add 1 for each direction patch -> at least one direction has to be chosen for this patch
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                    continue;

                if ( verbose && (lid1 == 0) ) std::cout << "[" << __func__ << "]: " << "Constraining " << prims[lid][lid1].getTag( _PrimitiveT::GID ) << " to choose one of ";
                coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = 1.0;
                non_zero_line = true;
                if ( verbose ) std::cout << prims[lid][lid1].getTag( _PrimitiveT::DIR_GID ) << ", ";
            }
            if ( verbose )std::cout << " directions";

            if ( non_zero_line )
            {
                // unique insertion
                unsigned prev_size = uniqueA.size(); // store prev set size to check for unique insert
                uniqueA.insert( coeffs );            // insert into set
                if ( uniqueA.size() > prev_size )    // if line was unique, add to problem
                {
                    problem.addLinConstraint( _OptProblemT::BOUND::GREATER_EQ, 1, problem.getINF(), coeffs ); // 1 <= A( lid, : ) * X <= INF
                    if  ( verbose )     std::cout << " ADDED\n";
                }
                else
                    if ( verbose )          std::cout << " IGNORED\n";
            }
        } // ... for each patch

        return err;
    } // ...everyPatchNeedsDirectionConstraint

    //! \brief                      Hybrid mode, large patches need one direction, points in small patches need to be assigned to one.
    //! \tparam _PrimitiveT         Concept: _PrimitiveContainerT::value_type::value_type
    //! \tparam _PointPrimitiveT    Concept: _PointContainerT::value_type
    //! \tparam _AssocT             Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
    //! \param[in,out] problem      The optimization problem to add to.
    template < class _PointPrimitiveDistanceFunctor
             , class _PrimitiveT
             , class _PointPrimitiveT
             , typename _Scalar
             , class _OptProblemT
             , class _PrimitiveContainerT
             , class _PointContainerT
             , class _AssocT
             , class _WeightsT
             >
    static inline int
    largePatchesNeedDirectionConstraint( _OptProblemT              & problem
                                       , _PrimitiveContainerT const& prims
                                       , _PointContainerT     const& points
                                       , _AssocT              const& lids_varids
                                       , _WeightsT            const& /*weights*/
                                       , _Scalar              const  scale
                                       , int                  const  pop_limit
                                       , bool                 const  verbose )
    {
        std::cerr << "[" << __func__ << "]: " << "don't use this" << std::endl;
        throw new std::runtime_error("don't use largePatchesNeedDirectionConstraint");

        typedef typename _AssocT::key_type      IntPair;        // <lid,lid1> pair to retrieve the linear index of a variable in the problem using lids_varids
        typedef typename _OptProblemT::Scalar   ProblemScalar;  // Scalar in OptProblem, usually has to be double because of the implementation
        typedef std::vector<ProblemScalar>      RowInA;         // "A" is the constraint matrix
        std::pair<int,int> stats(0,0); // debug

        int err = EXIT_SUCCESS;

        // count populations
        //typedef std::set<int>         PIdList;        // unique list of point ids
        //typedef std::map<int,PIdList> PopT;
        GidPidSetMap populations; // populations[gid] = list of points with GID==patch_id
        processing::getPopulations( populations, points );

        // Work on large patches, and record small patches
        // One direction / patch needs to be choosen for large patches
        RowInA             coeffs ( problem.getVarCount(), 0 ); // constraint coefficient line in A
        std::set< RowInA > uniqueA;                             // to ensure unique constraints
        PidSet             smalls_points;                       // contains point ids for points that are in small patches

        // for all patches
        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            if ( !prims[lid].size() ) continue;
            // cache patch group id
            const int gid = prims[lid].at(0).getTag(_PrimitiveT::GID);
            // check for population size
//            if ( populations[gid].size() < pop_limit )
//            { // small patch
//                if ( populations[gid].size() < 1 ) std::cerr << "[" << __func__ << "]: " << "patch with no points...needs to be expensive!! TODO" << std::endl;
//                // add all points to point-wise constraints todo list
//                PidSet::const_iterator end_it = populations[gid].end();
//                for ( PidSet::const_iterator it = populations[gid].begin(); it != end_it; ++it )
//                    smalls_points.insert( *it );
//            }
//            else
            if ( populations[gid].size() > pop_limit )
            { // create patch-wise constraint
                // init row to 0
                std::fill( coeffs.begin(), coeffs.end(), 0 );

                // add 1 for each direction patch -> at least one direction has to be chosen for this patch
                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                    coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = ProblemScalar( 1. );

                // unique insertion
                unsigned prev_size = uniqueA.size(); // store prev set size to check for unique insert
                uniqueA.insert( coeffs );            // insert into set
                if ( uniqueA.size() > prev_size )    // if line was unique, add to problem
                {
                    ++stats.first;
                    problem.addLinConstraint( _OptProblemT::BOUND::GREATER_EQ, 1, problem.getINF(), coeffs ); // 1 <= A( lid, : ) * X <= INF
                }
            }
        } // ... for each patch
//        // set-up point-wise constraints for small patches
//        PidSet::const_iterator end_it = smalls_points.end();
//        for ( PidSet::const_iterator it = smalls_points.begin(); it != end_it; ++it )
//        {
//            const int pid = *it;

//            // work - check distance to all lines, select the ones in range
//            RowInA coeffs  ( problem.getVarCount(), 0 );
//            bool   non_zero( false );
//            // for each patch and direction
//            for ( size_t lid = 0; lid != prims.size(); ++lid )
//                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
//                {
//                    _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );
//                    if ( dist < scale )
//                    {
//                        coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = ProblemScalar( 1. ) ;
//                        non_zero = true;
//                    }
//                }

//            if ( non_zero )
//            {
//                unsigned prev_size = uniqueA.size();
//                uniqueA.insert( coeffs );
//                if ( uniqueA.size() > prev_size )
//                {
//                    ++stats.second;
//                    problem.addLinConstraint( _OptProblemT::BOUND::GREATER_EQ, 1, problem.getINF(), coeffs );
//                }
//            }
//        } //... for each point
        if ( verbose ) std::cout << "[" << __func__ << "]: " << "added "  << stats.first << " large and " << stats.second << " small patch constraint lines" << std::endl;

        return err;
    } // ...everyPatchNeedsDirectionConstraint


    //____________________________________________DataCosts____________________________________________

    //! \brief              Adds unary costs to problem based on scale wide band assocation.
    //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
    template < class _PointPrimitiveDistanceFunctor
             , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
             , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
             , typename _Scalar
             , class _OptProblemT
             , class _PrimitiveContainerT
             , class _PointContainerT
             , class _AssocT
             , class _WeightsT
             >
    static inline int
    bandBasedDataCost( _OptProblemT              & problem
                     , _PrimitiveContainerT const& prims
                     , _PointContainerT     const& points
                     , _AssocT              const& lids_varids
                     , _WeightsT            const& weights
                     , _Scalar              const  scale )
    {
        typedef typename _AssocT::key_type IntPair;

        int err = EXIT_SUCCESS;

        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                    continue;

                // unary cost
                {
                    unsigned cnt = 0;
                    _Scalar unary_i = _Scalar(0);
                    for ( size_t pid = 0; pid != points.size(); ++pid )
                    {
                        // if within scale, add unary cost
                        _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );
                        if ( dist < scale ) //TODO: patch membership
                        {
                            unary_i += dist;
                            ++cnt;
                        }
                            //unary_i += dist / scale;
                    } // for points

                    _Scalar coeff = /* complx: */ weights(2) + /* unary: */ weights(0) * unary_i / _Scalar(cnt);
                    if ( !cnt )
                        coeff = weights(2) + weights(0) * 2;

                    problem.addLinObjective( /* var_id: */ lids_varids.at( IntPair(lid,lid1) )
                                           , coeff );
                }
            }
        } //...for patches

        return err;
    } //...bandBasedDataCost

    //! \brief Adds unary costs to problem based on point to primitive associations.
    //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
    template < class _PointPrimitiveDistanceFunctor
             , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
             , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
             , typename _Scalar
             , class _OptProblemT
             , class _PrimitiveContainerT
             , class _PointContainerT
             , class _AssocT
             , class _WeightsT
             >
    static inline int
    associationBasedDataCost( _OptProblemT              & problem
                            , _PrimitiveContainerT const& prims
                            , _PointContainerT     const& points
                            , _AssocT              const& lids_varids
                            , _WeightsT            const& weights
                            , _Scalar              const  /*scale*/
                            , _Scalar              const freq_weight
                            , bool                 const verbose )
    {
        typedef typename _AssocT::key_type IntPair;

        const _Scalar w_mod_base     = ProblemSetupParams<_Scalar>::w_mod_base,
                      w_mod_base_inv = _Scalar(1.) - w_mod_base;

        int err = EXIT_SUCCESS;

        // INSTANCES
        // added 17/09/2014 by Aron
        // modified 22/09/2014 by Aron
        std::map< int, int > dir_instances;
        int active_count = 0;
        for ( size_t lid = 0; lid != prims.size(); ++lid )
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) != _PrimitiveT::STATUS_VALUES::ACTIVE )
                    continue;

                ++dir_instances[ prims[lid][lid1].getTag(_PrimitiveT::DIR_GID) ];
                ++active_count;
            }

        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            // check, if any directions for patch
            if ( !prims[lid].size() )
            {
                //std::cerr << "[" << __func__ << "]: " << "no directions for patch[" << lid << "]! This shouldn't happen. Skipping patch..." << std::endl;
                continue;
            }

            // cache patch group id to match with point group ids
            const int gid = prims[lid][0].getTag( _PrimitiveT::GID );

            // for each direction
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                    continue;

                // point count for normalization
                unsigned cnt = 0;
                // data-cost coefficient (output)
                _Scalar unary_i = _Scalar(0);
                // for each point, check if assigned to main patch (TODO: move assignment test to earlier, it's indep of lid1)
                for ( size_t pid = 0; pid != points.size(); ++pid )
                {
                    //if ( points[pid].getTag( PointPrimitiveT::GID ) == static_cast<int>(lid) )
                    if ( points[pid].getTag( _PointPrimitiveT::GID ) == gid )
                    {
                        // if within scale, add unary cost
                        _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );
                        unary_i += dist * dist; //changed on 18/09/14
                        ++cnt;              // normalizer
                    }
                } // for points

                // average data cost
                _Scalar coeff = cnt ? /* unary: */ weights(0) * unary_i / _Scalar(cnt)
                                    : /* unary: */ weights(0) * _Scalar(2);            // add large weight, if no points assigned

                // prefer dominant directions
                if ( freq_weight > _Scalar(0.) )
                {
                    const int dir_gid = prims[lid][lid1].getTag( _PrimitiveT::GID );

                    if ( verbose && dir_instances[dir_gid] )
                        std::cout << "[" << __func__ << "]: " << "changed " << coeff << " to ";

                    // changed by Aron 10:32 24/09/2014
                    // old version: 1/#did, better version would be normalized, so: 1 / (#did/all)
                    // new version: Dataweight = dataweight * (.1 + .9 * ((#did/n)^2 - 1.)^6)
                    if ( dir_instances[dir_gid] > 0 )
                    {
                        _Scalar v = _Scalar(dir_instances[dir_gid]) / _Scalar(active_count); // #did/n
                        v *= v;                                                              // (#did/n)^2
                        v -= _Scalar(1.);                                                    // (#did/n)^2 - 1.
                        v *= v;                                                              // ((#did/n)^2 - 1.)^2
                        v *= v * v;                                                          // ((#did/n)^2 - 1.)^6
                        coeff *= (w_mod_base + w_mod_base_inv * v) * freq_weight;            // .1 + .9 * ((#did/n)^2 - 1.)^6
                        //coeff *= freq_weight * _Scalar(1.) / ( _Scalar(1.) + std::log(dir_instances[dir_gid]) );
                        //coeff *= freq_weight / _Scalar(dir_instances[dir_gid]);
                    }
                    else
                        coeff *= freq_weight;

                    if ( verbose && dir_instances[dir_gid] )
                        std::cout << coeff << " since dirpop: " << dir_instances[dir_gid] << std::endl;
                }

                // complexity cost:
                coeff += weights(2); // changed by Aron on 21/9/2014

                // add to problem
                problem.addLinObjective( /* var_id: */ lids_varids.at( IntPair(lid,lid1) )
                                       , /*  value: */ coeff );
            } //...for each direction
        } //...for each patch

        return err;
    } //...associationBasedDataCost

#if 0
    //! \brief Nic's version, unfinished!.
    //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
    //! \warning Unfinished
    template < class _PointPrimitiveDistanceFunctor
             , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
             , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
             , typename _Scalar
             , class _OptProblemT
             , class _PrimitiveContainerT
             , class _PointContainerT
             , class _AssocT
             , class _WeightsT
             >
    static inline int
    instanceBasedDataCost( _OptProblemT              & problem
                            , _PrimitiveContainerT const& prims
                            , _PointContainerT     const& points
                            , _AssocT              const& lids_varids
                            , _WeightsT            const& weights
                            , _Scalar              const  /*scale*/ )
    {
            int err = EXIT_SUCCESS;
            std::cerr << "[" << __func__ << "]: " << "UNIMPLEMENTED" << std::endl;
    //        std::map<int, std::pair<Scalar,int> > repr_dirs; // [lid] = <cost,npoints>
    //        // unary cost -> lin objective
    //        for ( size_t lid = 0; lid != prims.size(); ++lid )
    //        {
    //            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
    //            {
    //                unsigned cnt = 0;
    //                Scalar unary_i = Scalar(0);
    //                for ( size_t pid = 0; pid != points.size(); ++pid )
    //                {
    //                    if ( points[pid].getTag( PointPrimitiveT::GID ) == static_cast<int>(lid) )
    //                    {
    //                        // if within scale, add unary cost
    //                        Scalar dist = PointPrimitiveDistanceFunctor::eval<Scalar>( points[pid], prims[lid][lid1] );
    //                        unary_i += dist;
    //                        ++cnt;
    //                    }
    //                } // for points

    //                if ( prims[lid].getTag(PrimitiveT::GID) != lid )
    //                    std::cerr << "[" << __func__ << "]: " << "a;sldjfa;lkjdf, " << prims[lid].getTag(PrimitiveT::GID) << "lid != tag" << prims[lid].getTag(PrimitiveT::GID) << std::endl;

    //                if ( prims[lid].getTag(PrimitiveT::GID) == prims[lid1].getTag(PointPrimitiveT::GID) )
    //                    repr_dirs[ lid ] = std::pair<Scalar,int>( unary_i, cnt );
    //            }
    //        } // for lines

    //        for ( size_t lid = 0; lid != prims.size(); ++lid )
    //        {
    //            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
    //            {
    //                Scalar coeff = /* complx: */ weights(2) + /* unary: */ weights(0) * unary_i / Scalar(cnt);
    //                if ( !cnt )
    //                    coeff = weights(2) + weights(0) * 2;

    //                problem.addLinObjective( /* var_id: */ lids_varids[Key(lid,lid1)]
    //                        , coeff );
    //            }
    //        }
        return err;
    } //...instanceBasedDataCost()
#endif

} //...namespace ProblemSetup
} //...namespace GF2

#endif // GF2_PROBLEMSETUP_HPP
