#ifndef GF2_PROBLEMSETUP_HPP
#define GF2_PROBLEMSETUP_HPP

#include <iostream>                               // cout, cerr, endl
#include "Eigen/Dense"

#if GF2_USE_PCL
#   include "pcl/console/parse.h" // pcl::console::parse_argument
#endif

#include "qcqpcpp/optProblem.h"                   // OptProblem

#include "globfit2/optimization/energyFunctors.h" // AbstractPrimitivePrimitiveEnergyFunctor,
#include "globfit2/parameters.h"                  // ProblemSetupParams
#include "globfit2/processing/util.hpp"           // getPopulation()
#include "globfit2/io/io.h"                       // readPrimitives(), readPoints()

namespace GF2 {

//____________________________________________class ProblemSetup __________________________________________________

template < class _PrimitiveContainerT
         , class _PointContainerT
         , class _PrimitiveT
         , class _PointPrimitiveT
         > int
ProblemSetup::formulateCli( int    argc
                          , char** argv )
{
    bool verbose = false;
    typedef          MyPointPrimitiveDistanceFunctor               _PointPrimitiveDistanceFunctor;
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
            /**/ if ( !constr_mode_str.compare( "patch"  ) ) params.constr_mode = ProblemSetupParams<Scalar>::PATCH_WISE; // default
            else if ( !constr_mode_str.compare( "point"  ) ) params.constr_mode = ProblemSetupParams<Scalar>::POINT_WISE; // banded
            else if ( !constr_mode_str.compare( "hybrid" ) ) params.constr_mode = ProblemSetupParams<Scalar>::HYBRID; // hybrid
            else std::cerr << "[" << __func__ << "]: " << "could NOT parse constraint mode, assuming " << (int)params.constr_mode << std::endl;
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
                      << " [--cost-fn *" << cost_string << "* (cexp | sqrt)" << "]\n"
                      << " [--data-mode *" << (int)params.data_cost_mode << "* (assoc | band | instance) ]\n"
                      << " [--constr-mode *" << (int)params.constr_mode << "* (patch | point | hybrid ) ]\n"
                      << " [--srand " << srand_val << "]\n"
                      << " [--rod " << problem_rel_path << "]\t\tRelative output path of the output matrix files\n"
                      << " [--patch-pop-limit " << params.patch_population_limit << "]\n"
                      << " [--freq-weight" << params.freq_weight << "]\n"
                      << std::endl;
            if ( !verbose )
                return EXIT_FAILURE;
        }
    } // ... parse params

    // Read desired angles
    {
        processing::appendAnglesFromGenerators( params.angles, angle_gens, verbose );
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

    AbstractPrimitivePrimitiveEnergyFunctor<Scalar,_PrimitiveT> *primPrimDistFunctor = NULL;
    // parse cost function
    {
             if ( !cost_string.compare("cexp") )    primPrimDistFunctor = new CExpPrimitivePrimitiveEnergyFunctor<Scalar,_PrimitiveT>( params.angles );
        else if ( !cost_string.compare("sqrt") )    primPrimDistFunctor = new SqrtPrimitivePrimitiveEnergyFunctor<Scalar,_PrimitiveT>( params.angles );
        else                                        std::cerr << "[" << __func__ << "]: " << "Could not parse cost functor input..." << std::endl;
    } //...parse cost function

    // WORK
    problemSetup::OptProblemT problem;
    int err = formulate<_PointPrimitiveDistanceFunctor>( problem
                                                       , prims
                                                       , points
                                                       , params.constr_mode
                                                       , params.data_cost_mode
                                                       , params.scale
                                                       , params.weights
                                                       , primPrimDistFunctor
                                                       , params.patch_population_limit
                                                       , params.dir_id_bias
                                                       , !calc_energy && verbose
                                                       , params.freq_weight );

    // dump. default output: ./problem/*.csv; change by --rod
    if ( EXIT_SUCCESS == err )
    {
        if ( !calc_energy )
        {
            std::string parent_path = boost::filesystem::path(candidates_path).parent_path().string();
            if ( !parent_path.empty() )     parent_path += "/";
            else                            parent_path =  ".";
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
        }
    } //...dump

    // cleanup
    if ( primPrimDistFunctor ) { delete primPrimDistFunctor; primPrimDistFunctor = NULL; }

    return err;
} //...ProblemSetup::formulateCli()


template < class _PointPrimitiveDistanceFunctor
         , class _PrimitiveContainerT
         , class _PointContainerT
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
                       , AbstractPrimitivePrimitiveEnergyFunctor<_Scalar,_PrimitiveT>* const& primPrimDistFunctor
                       , int                                                           const  patch_pop_limit
                       , _Scalar                                                       const  dir_id_bias
                       , int                                                           const  verbose
                       , _Scalar                                                       const  freq_weight /* = 0. */
        )
{
    using problemSetup::OptProblemT;

    // work - formulate problem
    int err = EXIT_SUCCESS;

    // log
    if ( verbose ) { std::cout << "[" << __func__ << "]: " << "formulating problem...\n"; fflush(stdout); }

    typedef std::pair<int,int> IntPair;
    /*   */ std::map <IntPair,int> lids_varids;
    std::set<int> chosen_varids;

    // ____________________________________________________
    // Variables - Add variables to problem
    {
        char name[16];
        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                // add var
                int gid     = prims[lid][lid1].getTag( _PrimitiveT::GID );
                int dir_gid = prims[lid][lid1].getTag( _PrimitiveT::DIR_GID );
                sprintf( name, "x_%d_%d", gid, dir_gid );

                // store var_id for later, add binary variable
                const int var_id = problem.addVariable( OptProblemT::BOUND::RANGE, 0.0, 1.0, OptProblemT::VAR_TYPE::INTEGER );
                lids_varids[ IntPair(lid,lid1) ] = var_id;
                if ( prims[lid][lid1].getTag(_PrimitiveT::CHOSEN) > 0 )
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
                        ( problem, prims, points, lids_varids, weights, scale );
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
                        ( problem, prims, points, lids_varids, weights, scale, freq_weight );
                break;

            case ProblemSetupParams<_Scalar>::DATA_COST_MODE::INSTANCE_BASED:
                err = problemSetup::instanceBasedDataCost<_PointPrimitiveDistanceFunctor, _PrimitiveT, _PointPrimitiveT>
                        ( problem, prims, points, lids_varids, weights, scale );
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
        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                for ( size_t olid = 0; olid != prims.size(); ++olid )
                {
                    for ( size_t olid1 = 0; olid1 != prims[olid].size(); ++olid1 )
                    {
                        // skip same line, that's always zero
                        if ( (lid == olid) && (lid1 == olid1) ) continue;

                        _Scalar dist = primPrimDistFunctor->eval( prims[lid][lid1], prims[olid][olid1] );
                        if ( prims[lid][lid1].getTag(_PrimitiveT::DIR_GID) != prims[olid][olid1].getTag(_PrimitiveT::DIR_GID) )
                            dist += dir_id_bias;
                        dist *= weights(1);
                        problem.addQObjective( lids_varids[ IntPair(lid,lid1) ], lids_varids[ IntPair(olid,olid1) ], dist );
                    } // ... olid1
                } // ... olid
            } // ... lid1
        } // ... lid
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
                                      , _Scalar              const  /*scale*/ )
    {
        bool verbose = false;
        typedef typename _AssocT::key_type IntPair;

        int err = EXIT_SUCCESS;

        // one direction / patch needs to be choosen
        std::vector< double              > coeffs ( problem.getVarCount(), 0 ); // constraint coefficient line in A
        std::set   < std::vector<double> > uniqueA;                             // to ensure unique constraints
        // for all patches
        for ( size_t lid = 0; lid != prims.size(); ++lid )
        {
            // init row to 0
            std::fill( coeffs.begin(), coeffs.end(), 0 );

            // add 1 for each direction patch -> at least one direction has to be chosen for this patch
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                if ( verbose && (lid1 == 0) ) std::cout << "[" << __func__ << "]: " << "Constraining " << prims[lid][lid1].getTag( _PrimitiveT::GID ) << " to choose one of ";
                coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = 1.0;
                if ( verbose ) std::cout << prims[lid][lid1].getTag( _PrimitiveT::DIR_GID ) << ", ";
            }
            if ( verbose )std::cout << " directions";

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
            if ( populations[gid].size() < pop_limit )
            { // small patch
                if ( populations[gid].size() < 1 ) std::cerr << "[" << __func__ << "]: " << "patch with no points...needs to be expensive!! TODO" << std::endl;
                // add all points to point-wise constraints todo list
                PidSet::const_iterator end_it = populations[gid].end();
                for ( PidSet::const_iterator it = populations[gid].begin(); it != end_it; ++it )
                    smalls_points.insert( *it );
            }
            else
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

        // set-up point-wise constraints for small patches
        PidSet::const_iterator end_it = smalls_points.end();
        for ( PidSet::const_iterator it = smalls_points.begin(); it != end_it; ++it )
        {
            const int pid = *it;

            // work - check distance to all lines, select the ones in range
            RowInA coeffs  ( problem.getVarCount(), 0 );
            bool   non_zero( false );
            // for each patch and direction
            for ( size_t lid = 0; lid != prims.size(); ++lid )
                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                {
                    _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );
                    if ( dist < scale )
                    {
                        coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = ProblemScalar( 1. ) ;
                        non_zero = true;
                    }
                }

            if ( non_zero )
            {
                unsigned prev_size = uniqueA.size();
                uniqueA.insert( coeffs );
                if ( uniqueA.size() > prev_size )
                {
                    ++stats.second;
                    problem.addLinConstraint( _OptProblemT::BOUND::GREATER_EQ, 1, problem.getINF(), coeffs );
                }
            }
        } //... for each point

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
                            , _Scalar              const freq_weight )
    {
        typedef typename _AssocT::key_type IntPair;

        int err = EXIT_SUCCESS;

        // INSTANCES // added 17/09/2014 by Aron
        std::map< int, int > dir_instances;
        for ( size_t lid = 0; lid != prims.size(); ++lid )
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                ++dir_instances[ prims[lid][lid1].getTag(_PrimitiveT::DIR_GID) ];
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
                        unary_i += dist;
                        ++cnt;              // normalizer
                    }
                } // for points

                _Scalar coeff = cnt ? /* complx: */ weights(2) + /* unary: */ weights(0) * unary_i / _Scalar(cnt)
                                    : /* complx: */ weights(2) + /* unary: */ weights(0) * _Scalar(2);            // add large weight, if no points assigned
                if ( freq_weight > _Scalar(0.) )
                {
                    const int dir_gid = prims[lid][lid1].getTag( _PrimitiveT::GID );

                    std::cout << "[" << __func__ << "]: " << "changed " << coeff << " to ";
                    if ( dir_instances[dir_gid] > 0 )
                        coeff *= freq_weight * _Scalar(1.) / _Scalar(dir_instances[dir_gid]);
                    std::cout << coeff << " since dirpop: " << dir_instances[dir_gid] << std::endl;
                }

                // add to problem
                problem.addLinObjective( /* var_id: */ lids_varids.at( IntPair(lid,lid1) )
                                       , /*  value: */ coeff );
            } //...for each direction
        } //...for each patch

        return err;
    } //...associationBasedDataCost

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

} //...namespace ProblemSetup
} //...namespace GF2

#endif // GF2_PROBLEMSETUP_HPP
