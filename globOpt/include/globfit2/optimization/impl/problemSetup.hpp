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
#include "globfit2/processing/angle_util.hpp"     // appendAngle...
#include "globfit2/io/io.h"                       // readPrimitives(), readPoints()

#include "globfit2/processing/graph.hpp"
#include "globfit2/processing/angle_util.hpp" // appendAnglefromgen
#include "omp.h"

namespace GF2 {

//____________________________________________class ProblemSetup __________________________________________________

template <typename _Scalar, class _PrimitiveT, class _AnglesT>
inline _Scalar calcPwCost( _PrimitiveT const& p0, _PrimitiveT const& p1, _AnglesT const& angles )
{
    //weights(1) * sqrt( MyPrimitivePrimitiveAngleFunctor::eval( *p0, *p1, angles ) ); //changed 16:11 15/01/2015
    return std::sqrt( MyPrimitivePrimitiveAngleFunctor::eval( p0, p1, angles ) );
}

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
    AnglesT                   angle_gens( {AnglesT::Scalar(90.)} );
    int                       srand_val          = 123456;
    std::string               cost_string        = "sqrt";
    std::string               problem_rel_path   = "problem";
    std::string               data_cost_mode_str = "assoc";
    std::string               constr_mode_str    = "hybrid";
    std::string               energy_path        = "energy.csv";
    int                       clustersMode       = 1;
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
        Scalar collapseAngleDeg = params.collapseAngleSqrt * params.collapseAngleSqrt * Scalar( 180. ) / M_PI;
        pcl::console::parse_argument( argc, argv, "--collapse-angle-deg", collapseAngleDeg );
        params.collapseAngleSqrt = sqrt( collapseAngleDeg * M_PI / Scalar(180.) );

        // energy
        {
            calc_energy = pcl::console::find_switch( argc, argv, "--energy" );
        }

        // freq_weight
        pcl::console::parse_argument( argc, argv, "--freq-weight", params.freq_weight );

        // clusters
        clustersMode = !pcl::console::find_argument( argc, argv, "--no-clusters" );

        pcl::console::parse_argument( argc, argv, "--spat-weight", params.spatial_weight_coeff );
        pcl::console::parse_argument( argc, argv, "--spat-dist-mult", params.spatial_weight_dist_mult );

        params.useAngleGen = pcl::console::find_switch( argc, argv, "--use-angle-gen" );
        pcl::console::parse_argument( argc, argv, "--trunc-angle", params.truncAngle );

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
                      << " [--freq-weight " << params.freq_weight << "]\n"
                      << " [--energy-out " << energy_path << "]\n"
                      << " [--no-paral]\n"
                      << " [--no-clusters " << clustersMode << "]\n"
                      << " [--spat-weight " << params.spatial_weight_coeff << "]\t How much penalty is added for mismatching primitives for patches in proximity.\n"
                      << " [--spat-dist-mult" << params.spatial_weight_dist_mult << "]\t How many times scale is proximity threshold \n"
                      << " [--use-angle-gen " << params.useAngleGen << "]\n"
                      << " [--trunc-angle " << params.truncAngle << "]\n"
                      << " [--collapse-angle-deg " << ((params.collapseAngleSqrt*params.collapseAngleSqrt) * 180. / M_PI)  << "]\n"
                      << std::endl;
            if ( !verbose )
                return EXIT_FAILURE;
        }
    } // ... parse params

    // Read desired angles
    bool no_paral = pcl::console::find_switch( argc, argv, "--no-paral");
    {
        angles::appendAnglesFromGenerators( params.angles, angle_gens, no_paral, verbose );
    } //...read angles

    // read points
    _PointContainerT     points;

    PclCloudPtrT pclCloud( new PclCloudT() );
    {
        if ( verbose ) std::cout << "[" << __func__ << "]: " << "reading cloud from " << cloud_path << "...";
        io::readPoints<_PointPrimitiveT>( points, cloud_path, &pclCloud );
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
    std::vector<std::pair<PidT,LidT> > points_primitives;
    io::readAssociations( points_primitives, assoc_path, NULL );
    int nrWarnings = 0;
    for ( size_t i = 0; i != points.size(); ++i )
    {
        if ( i > points_primitives.size() )
        {
            std::cerr << "more points than associations..." << std::endl;
            return EXIT_FAILURE;
        }

        points[i].setTag( _PointPrimitiveT::TAGS::GID, points_primitives[i].first );

        if ( points[i].getTag(_PointPrimitiveT::TAGS::GID) == -1 )
            ++nrWarnings;
    } //...read associations
    if ( nrWarnings )
        std::cerr << "[" << __func__ << "]: " << nrWarnings << "/" << points.size() << " points unAssigned!" << std::endl;

    //AbstractPrimitivePrimitiveEnergyFunctor<Scalar,_PrimitiveT> *primPrimDistFunctor = NULL;
    SpatialSqrtPrimitivePrimitiveEnergyFunctor<_FiniteFiniteDistFunctor, _PointContainerT, Scalar,_PrimitiveT> *primPrimDistFunctor = NULL;
    // parse cost function
    {
        if ( !cost_string.compare("spatsqrt") )
        {
            primPrimDistFunctor = new SpatialSqrtPrimitivePrimitiveEnergyFunctor<_FiniteFiniteDistFunctor, _PointContainerT, Scalar,_PrimitiveT>
                    ( params.angles, points, params.scale );
            primPrimDistFunctor->_verbose = verbose;
            primPrimDistFunctor->setUseAngleGen( params.useAngleGen );
            primPrimDistFunctor->setDirIdBias  ( params.dir_id_bias );
            primPrimDistFunctor->setTruncAngle ( params.truncAngle );
            primPrimDistFunctor->setSpatialWeightCoeff( params.spatial_weight_coeff );
            primPrimDistFunctor->setSpatialWeightDistMult( params.spatial_weight_dist_mult );
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
#if 1
    int err = formulate2<_PointPrimitiveDistanceFunctor>( problem
                                                        , prims
                                                        , points
                                                        , params.constr_mode
                                                        , params.data_cost_mode
                                                        , params.scale
                                                        , params.weights
                                                        , primPrimDistFunctor
                                                        , angle_gens_in_rad
                                                        , params.patch_population_limit
                                                        , pclCloud
                                                        , !calc_energy && verbose
                                                        , params.freq_weight
                                                        , clustersMode
                                                        , params.collapseAngleSqrt
                                                        );
#else
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
                                                           //, params.dir_id_bias
                                                           , !calc_energy && verbose
                                                           , params.freq_weight
                                                           , clustersMode );
#endif

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

/*! \brief                  Calculate vicinity of patches based on smallest point-point distance.
 * \tparam      NeighMapT   map<GidT,set<GidT>>
 * \param[in]   radius      Lookup radius, usually 2x scale (\ref ProblemSetupParams::spatial_weight_distance)
 */
template <class NeighMapT, typename _PointContainerT, typename _Scalar>
inline void calculateNeighbourhoods( NeighMapT &proximity, _PointContainerT const& points, const _Scalar radius )
{
    //typedef typename _PointContainerT::value_type PointPrimitiveT;
    //typedef typename PointPrimitiveT::Scalar Scalar;
    using pclutil::PclSearchPointT;
    typedef Eigen::Vector3f Colour;

    pclutil::PclSearchTreePtrT tree = pclutil::buildANN( points );

    std::vector<int>    k_indices;
    std::vector<float>  k_sqr_distances;
    int warningCount = 0;
    #pragma omp parallel for private(k_indices,k_sqr_distances) num_threads(GF2_MAX_OMP_THREADS)
    for ( size_t i = 0; i < points.size(); ++i )
    {
        const GidT gidI = points[i].getTag(PointPrimitiveT::TAGS::GID);

        if ( gidI == PointPrimitiveT::TAG_UNSET ) continue;

        pclutil::PclSearchPointT pnt;
        pnt.getVector3fMap() = points[i].template pos();
        tree->radiusSearch( pnt, radius, k_indices, k_sqr_distances, /*maxnn:*/ 0 );

#       pragma omp critical (NEIGHWARN)
        if ( k_indices.size() > 1000 )
        {
            ++warningCount;
        }

        for ( size_t j = 1; j < k_indices.size(); ++j )
        {
            const PidT neighPid = k_indices[j];
            const GidT gidJ = points[neighPid].getTag(PointPrimitiveT::TAGS::GID);
            if (    ( gidJ == PointPrimitiveT::TAG_UNSET )
                 || ( gidJ == gidI )
               ) continue;

#           pragma omp critical (PROXIMITY)
            {
                proximity[ gidI ].insert( gidJ );
                proximity[ gidJ ].insert( gidI );
            }
        } //...foreach neighbour

    } //...foreach point
    std::cerr << "[" << __func__ << "]: " << "more, than 1000 neighbrours " << warningCount << "/" << points.size() << " times" << std::endl;
} //...calculateNeighbourhoods

template < class _PointPrimitiveDistanceFunctor
         , class _PrimitiveContainerT
         , class _PointContainerT
         , class _PrimPrimDistFunctorT
         , class _PrimitiveT
         , class _PointPrimitiveT
         , typename _Scalar
         > int
ProblemSetup::formulate2( problemSetup::OptProblemT                                          & problem
                       , _PrimitiveContainerT                                          const& prims
                       , _PointContainerT                                              const& points
                       , typename ProblemSetupParams<_Scalar>::CONSTR_MODE             const  constr_mode
                       , typename ProblemSetupParams<_Scalar>::DATA_COST_MODE          const  data_cost_mode
                       , _Scalar                                                       const  scale
                       , Eigen::Matrix<_Scalar,-1,1>                                   const& weights
                       , _PrimPrimDistFunctorT                                       * const& primPrimDistFunctor // AbstractPrimitivePrimitiveEnergyFunctor<_Scalar,_PrimitiveT>
                       , AnglesT                                                       const& angle_gens_in_rad
                       , int                                                           const  patch_pop_limit
                       , PclCloudPtrT                                                            & pclCloud
                       , int                                                           const  verbose
                       , _Scalar                                                       const  freq_weight /* = 0. */
                       , int                                                           const  clusterMode
                       , _Scalar                                                       const  collapseThreshold /* = 0.07 */ // sqrt( 0.1 * PI / 180 ) == 0.06605545496
        )
{
    using problemSetup::OptProblemT;

    typedef Graph< _Scalar, typename MyGraphConfig<_Scalar>::UndirectedGraph > GraphT;
    typedef graph::EdgeT<_Scalar> EdgeT;
    typedef typename OptProblemT::SparseMatrix        SparseMatrix;
    typedef typename OptProblemT::SparseEntry         SparseEntry;

    bool needPairwise = ((primPrimDistFunctor->getSpatialWeightCoeff() != _Scalar(0.)) || clusterMode);
    if ( needPairwise ) std::cout << "needPairwise, because spatW: " << primPrimDistFunctor->getSpatialWeightCoeff() << ", and clusterMode: " << clusterMode << std::endl;

    // work - formulate problem
    int err = EXIT_SUCCESS;

    // log
    if ( verbose ) { std::cout << "[" << __func__ << "]: " << "formulating problem...\n"; fflush(stdout); }

    typedef std::pair<LidT   ,LidT> IntPair;
    /*   */ std::map <IntPair,LidT> lids_varids;
    /*   */ std::map <DidT   ,LidT> dIdsVarIds;
    std::set<LidT> chosen_varids;

    AnglesT angles = primPrimDistFunctor->getAngles();
    if ( angle_gens_in_rad.end() != std::find_if( angle_gens_in_rad.begin(), angle_gens_in_rad.end(), [](Scalar const& angle) { return angle > 2. * M_PI; } ) )
    {
        std::cerr << "[" << __func__ << "]: " << "angle_gens need to be in rad, are you sure?" << std::endl;
        std::cerr<<"angle_gens_in_rad:";for(size_t vi=0;vi!=angle_gens_in_rad.size();++vi)std::cerr<<angle_gens_in_rad[vi]<<" ";std::cerr << "\n";
        throw new std::runtime_error("angle_gens need to be in rad, are you sure");
    }

    GidPidVectorMap populations;
    processing::getPopulations( populations, points );

    // find smallest pwcost
    std::cout << "[" << __func__ << "]: " << "collapse loop start..." << std::endl; fflush(stdout);
    typedef std::pair<DidT,DidT> DIdPair;
    DIdPair minPair;
    _Scalar minScore = std::numeric_limits<_Scalar>::max();
    std::map< DidT, ULidT > dIdPopuls; // <did, pointcount>
    {
        std::set< DIdPair > visited;
        for ( size_t lid = 0; lid != prims.size(); ++lid )
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                // skip small
                //if ( prims[lid][lid1].getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL ) continue;
                if ( !(   (prims[lid][lid1].getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::ACTIVE)
                       || (prims[lid][lid1].getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::FIXED)
                      )
                   )
                    continue;

                const DidT did = prims[lid][lid1].getTag(_PrimitiveT::TAGS::DIR_GID);
                const GidT gid = prims[lid][lid1].getTag(_PrimitiveT::TAGS::GID);

                // skip empty
                if ( (populations.find(gid) != populations.end()) && !populations[gid].size() ) continue;

                dIdPopuls[did] += ( (populations.find( gid ) != populations.end()) ? populations[gid].size()
                                                                                   : 0 );

                for ( size_t lid2 = 0; lid2 != prims.size(); ++lid2 )
                    for ( size_t lid3 = 0; lid3 != prims[lid2].size(); ++lid3 )
                    {
                        //if ( prims[lid2][lid3].getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL ) continue;

                        if ( !(   (prims[lid2][lid3].getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::ACTIVE)
                               || (prims[lid2][lid3].getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::FIXED)
                              )
                           )
                            continue;

                        const DidT did2 = prims[lid2][lid3].getTag(_PrimitiveT::TAGS::DIR_GID);
                        const GidT gid2 = prims[lid2][lid3].getTag(_PrimitiveT::TAGS::GID);
                        if ( did == did2 ) continue;
                        if ( (populations.find(gid2) != populations.end()) && !populations[gid2].size() ) continue;

                        DIdPair pair = DIdPair(did,did2);
                        if ( visited.find( pair ) == visited.end() )
                        {
                            _Scalar score = calcPwCost<_Scalar>( prims[lid][lid1], prims[lid2][lid3], angles );
                            if ( score < minScore )
                            {
                                minScore = score;
                                minPair = pair;
                            }
                            visited.insert( pair );
                        } //...if not visited
                    } //...inner
            } //...outer
    } //...smallest pwcost
    std::cout << "[" << __func__ << "]: " << "mincost: " << minScore << " by " << minPair.first << "-" << minPair.second
              << ", populs: " << dIdPopuls[ minPair.first ] << " vs " << dIdPopuls[ minPair.second ]
              << std::endl;
    DIdPair replaceBy( _PrimitiveT::LONG_VALUES::UNSET, _PrimitiveT::LONG_VALUES::UNSET ); // replace first by second
    if ( minScore < collapseThreshold ) // sqrt( 0.1 * PI / 180 ) == 0.06605545496
    {
        replaceBy = minPair;
        if ( dIdPopuls[minPair.first] > dIdPopuls[minPair.second] )
            std::swap( replaceBy.first, replaceBy.second );
        if (  dIdPopuls[replaceBy.first] == 0 )
            replaceBy.first = replaceBy.second = _PrimitiveT::LONG_VALUES::UNSET;
        else
            std::cout << "[" << __func__ << "]: " << "replace " << replaceBy.first << " by " << replaceBy.second << std::endl;
    }
    else
        std::cout << "[" << __func__ << "]: " << "not replacing, since minscore: " << minScore << "(" << minPair.first << ","  << minPair.second << ")" << std::endl;
    std::cout << "[" << __func__ << "]: " << "collapse loop finish..." << std::endl; fflush(stdout);


    // ____________________________________________________
    // Variables - Add variables to problem
    //std::map< DidT, LidT > dIdsVarIds; // holds the cluster variable Ids for each direction Id

    std::cout << "[" << __func__ << "]: " << "var loop start..." << std::endl; fflush(stdout);
    /* { < didT, [ varId, varId, ...] >, ... }
     * Holds the variable Ids of primitives that will have to
     * be connected to an extra node representing this dId
     */
    std::map< DidT, std::vector< LidT > > dIdsPrimVarIds;
    std::map< DidT, _PrimitiveT const* > dIdsPrims; // representatives for direction read later
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
                GidT gId = prims[lid][lid1].getTag( _PrimitiveT::TAGS::GID     );
                LidT dId = prims[lid][lid1].getTag( _PrimitiveT::TAGS::DIR_GID );
                sprintf( name, "x_%d_%d", gId, dId );

                // store var_id for later, add binary variable
                const UidT var_id = problem.addVariable( OptProblemT::BOUND::RANGE, 0.0, 1.0, OptProblemT::VAR_TYPE::INTEGER
                                                      , OptProblemT::LINEARITY::LINEAR, name ); // changed to nonlinear by Aron on 29.12.2014
                lids_varids[ IntPair(lid,lid1) ] = var_id;

                // save for initial starting point: 1. [active && ( no replace OR has not dId to be replaced )] OR [ replace set && has the dId to replace by ]
                if (    (    (prims[lid][lid1].getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::ACTIVE)
                          &&
                             ( (replaceBy.first == _PrimitiveT::LONG_VALUES::UNSET) || (replaceBy.first != dId)       )
                        )
                     ||
                        ( (replaceBy.first != _PrimitiveT::LONG_VALUES::UNSET) && (replaceBy.second == dId)           )
                   )
                {
                    chosen_varids.insert( var_id );
                }

                // save dId node
                dIdsPrimVarIds[ dId ].push_back( var_id );
                //std::cout<<"dIdsPrimVarIds["<<dId<<"]:";for(size_t vi=0;vi!=dIdsPrimVarIds[dId].size();++vi)std::cout<<dIdsPrimVarIds[dId][vi]<<" ";std::cout << "\n";

                // save direction
                if ( dIdsPrims.find(dId) == dIdsPrims.end() )
                    dIdsPrims[ dId ] = &(prims[lid][lid1]);
            } //...for inner prims
        } //...for outer prims

        // add all dId nodes to the problem, and register them
        for ( auto dIdsIt = dIdsPrimVarIds.begin(); dIdsIt != dIdsPrimVarIds.end(); ++dIdsIt )
        {
            // first: dId
            // second: std::vector< varId >
            DidT dId     = dIdsIt->first;
            auto &varIds = dIdsIt->second;

            sprintf( name, "dId%uL", dId );

            // store varId for later, add binary variable
            const LidT varId = problem.addVariable( OptProblemT::BOUND::RANGE, 0.0, 1.0, OptProblemT::VAR_TYPE::BINARY
                                                 , OptProblemT::LINEARITY::NON_LINEAR, name );
            // record
            dIdsVarIds[ dId ] = varId;

            // select as initial, if any of it's primitives were selected
            for ( auto varIt = varIds.begin(); varIt != varIds.end(); ++varIt )
                if ( chosen_varids.find(*varIt) != chosen_varids.end() )
                {
                    chosen_varids.insert( varId );
                    break;
                }
        } //...for dIds
    } // ... variables
    std::cout << "[" << __func__ << "]: " << "var loop end..." << std::endl; fflush(stdout);

    std::cout << "[" << __func__ << "]: " << "constr loop start..." << std::endl; fflush(stdout);
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
                std::cerr << "[" << __func__ << "]: " << "POINT_WISE not used" << std::endl;
                throw new std::runtime_error("POINT_WISE not used");
                break;
            case ProblemSetupParams<_Scalar>::CONSTR_MODE::HYBRID:
                std::cerr << "[" << __func__ << "]: " << "HYBRID not used" << std::endl;
                throw new std::runtime_error("HYBRID not used");
                break;
        } //...switch constr_mode

        // Error check - Lin constraints
        if ( EXIT_SUCCESS != err )
        {
            std::cerr << "[" << __func__ << "]: " << "lin constraints setup returned error " << err << std::endl;
            return err;
        }
    } //...Lin constraints
    std::cout << "[" << __func__ << "]: " << "constr loop end..." << std::endl; fflush(stdout);
    std::cout << "[" << __func__ << "]: " << "data loop start..." << std::endl; fflush(stdout);
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
//                err = problemSetup::bandBasedDataCost<_PointPrimitiveDistanceFunctor, _PrimitiveT, _PointPrimitiveT>
//                        ( problem, prims, points, lids_varids, weights, scale );
                std::cerr << "::BAND_BASED not implemented" << std::endl;
                throw new std::runtime_error("::BAND_BASED not implemented");
                break;
        } //...switch data_cost_mode
    } //...unary cost

    // Error check - Unary costs
    if ( EXIT_SUCCESS != err )
    {
        std::cerr << "[" << __func__ << "]: " << "data cost setup returned error " << err << std::endl;
        return err;
    }
    std::cout << "[" << __func__ << "]: " << "data loop end..." << std::endl; fflush(stdout);

    // ____________________________________________________
    // Pairwise cost -> quad objective
    if ( EXIT_SUCCESS == err )
    {
        typedef std::vector< Eigen::Matrix<_Scalar,3,1> > ExtremaT;
        typedef std::pair<LidT,LidT>                      LidLid;
        typedef std::map< LidLid, ExtremaT >              ExtremaMapT;
        typedef OptProblemT::Scalar                       ProblemScalar;
        typedef typename GraphT::ComponentSizesT          ComponentSizesT;
        typedef typename GraphT::ComponentListT           ComponentListT;

        typedef std::map<GidT, std::set<GidT> >           ProximityMapT;
        ProximityMapT proximities;
        if ( needPairwise )
        {
            std::cout << "[" << __func__ << "]: " << "proximity start..." << std::endl; fflush(stdout);
            calculateNeighbourhoods( proximities, points, primPrimDistFunctor->getSpatialWeightDistMult() * scale );
            std::cout << "[" << __func__ << "]: " << "proximity end..." << std::endl; fflush(stdout);
        }

        //GraphT::testGraph();
        std::set< EdgeT > edgesList;

        const _Scalar halfSpatialWeightCoeff = primPrimDistFunctor->getSpatialWeightCoeff() / _Scalar(2.); // we add both ways, so adds up
        if ( needPairwise )
        {
            std::cout << "[" << __func__ << "]: " << "spatial start..." << std::endl; fflush(stdout);

#           pragma omp parallel for num_threads(GF2_MAX_OMP_THREADS)
            for ( size_t lid = 0; lid < prims.size(); ++lid )
            {
                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                {
                    _PrimitiveT const& prim = prims[lid][lid1];
                    if ( prim.getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                        continue;

                    const GidT gid = prim.getTag( _PrimitiveT::TAGS::GID );
                    const DidT did = prim.getTag( _PrimitiveT::TAGS::DIR_GID );

                    // extremas key
                    LidLid lidLid1( lid, lid1 );

                    for ( size_t lidOth = 0; lidOth != prims.size(); ++lidOth )
                    {
                        for ( size_t lid1Oth = 0; lid1Oth != prims[lidOth].size(); ++lid1Oth )
                        {
                            _PrimitiveT const& prim1 = prims[lidOth][lid1Oth];

                            if ( prim1.getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                                continue;

                            // skip same line, that's always zero
                            if ( (lid == lidOth) && (lid1 == lid1Oth) ) continue;

                            const GidT gIdOther = prim1.getTag( _PrimitiveT::TAGS::GID );
                            const DidT dIdOther = prim1.getTag( _PrimitiveT::TAGS::DIR_GID );

                            ProximityMapT::const_iterator gidNeighsIt = proximities.find( gid );
                            if (    ( did != dIdOther )
                                    && ( gid != gIdOther ) // we don't want to pollute problem with unnecessary edges
                                    && (    (    (gidNeighsIt != proximities.end()                               )
                                                 && (gidNeighsIt->second.find(gIdOther) != gidNeighsIt->second.end()) )
                                            //                                      || ()
                                            )
                                    )
                            {
                                //std::cout << "adding spatw " << halfSpatialWeightCoeff << " to " << gid << "-" << gIdOther << std::endl;
                                const LidT varId0 = lids_varids.at( lidLid1 );
                                const LidT varId1 = lids_varids.at( IntPair(lidOth,lid1Oth) );
#                               pragma omp critical (PS_PROBLEM)
                                {
                                    problem.addQObjective( varId0, varId1, halfSpatialWeightCoeff ); // /2, since it's going to be added both ways Aron 6/1/2015
                                }
                            }
                        } // ... olid1
                    } // ... olid
                } // ... lid1
            } // ... lid

            if ( clusterMode )
            {
                throw new std::runtime_error("turn off clusterMode!");

                GraphT graph( lids_varids.size() );
                for ( auto it = edgesList.begin(); it != edgesList.end(); ++it )
                {
                    graph.addEdge( it->_v0, it->_v1, /* not used right now: */ it->_w );
                }

                for ( LidT i = 0; i != lids_varids.size(); ++i )
                {
                    if ( !problem.getVarName(i).empty() )
                        graph.addVertexName( i, problem.getVarName(i) );
                }

                {
                    std::ofstream f;
                    f.open( "components.gv" );
                    f << "graph {\n";
                    ComponentListT components;
                    ComponentSizesT compSizes;
                    graph.getComponents( components, &compSizes );

                    std::map< PidT, std::vector<PidT> > clusters; // [ cluster0: [v0, v10,...], cluster1: [v3, v5, ...], ... ]
                    for ( LidT varId = 0; varId != components.size(); ++varId )
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

                    if ( clusterMode && clusters.size() )
                    {
                        //system( "dot -Tpng -o comps.png components.gv && (eog comps.png &)" );

                        // work
                        typename OptProblemT::SparseMatrix cluster_constraint( 1, problem.getVarCount()+1 );
                        char name[255];
                        LidT clusterId = 0;
                        for ( auto it = clusters.begin(); it != clusters.end(); ++it, ++clusterId )
                        {
                            sprintf(name,"cl_%d", clusterId );
                            LidT varid = problem.addVariable( OptProblemT::BOUND::RANGE, 0.0, 1.0, OptProblemT::VAR_TYPE::INTEGER
                                                             , OptProblemT::LINEARITY::LINEAR, name );

                            for ( LidT i = 0; i != it->second.size(); ++i )
                            {
                                cluster_constraint.insert( 0, it->second[i] ) = -1;
                            }
                            cluster_constraint.insert( 0, varid ) = (int)it->second.size();
                            // k * X_cluster_l <= A( l, : ) * X <= INF
                            problem.addConstraint( /*        type: */ OptProblemT::BOUND::GREATER_EQ
                                                 , /* lower_limit: */ 0 // k * X_cluster_l
                                                 , /* upper_limit: */ problem.getINF()
                                                 , /*      coeffs: */ &cluster_constraint );

                            chosen_varids.insert( varid );

                            std::cout << "added " << name << std::endl;

                            // extra!
                            //problem.addLinObjective( varid, -100 );
                        }
                    }
                }

                graph.draw( "graph.gv" );
                //system( "dot -Tpng -o graph.png graph.gv && (eog graph.png &)" );
            } //...if clusterMode
            std::cout << "[" << __func__ << "]: " << "spatial end..." << std::endl; fflush(stdout);
        } //...if spatialweight or clustermode
    } //...pairwise cost


    std::cout << "[" << __func__ << "]: " << "lvl2 constr start..." << std::endl; fflush(stdout);
    // ____________________________________________________
    // New dId constraint
    {
        for ( auto dIdIt = dIdsPrimVarIds.begin(); dIdIt != dIdsPrimVarIds.end(); ++dIdIt )
        {
            const LidT constraintId = problem.getConstraintCount();
            // first: dId
            //        didsVarIds[didIt->first] is the varId for this direction
            // second: vector< varId > of primitives with that Id

            SparseMatrix cluster_constraint( 1, problem.getVarCount() );
            std::vector<SparseEntry> quadConstraints;

            // add 1 for each primitive variable
            for ( LidT i = 0; i != dIdIt->second.size(); ++i )
            {
                // dIdIt->second[i]: prim_i

                cluster_constraint.insert( 0, dIdIt->second[i] ) = -1;
                quadConstraints.push_back( SparseEntry(dIdIt->second[i], dIdsVarIds[dIdIt->first], 1.) );
            }

            problem.addConstraint( /*        type: */ OptProblemT::BOUND::GREATER_EQ
                                 , /* lower_limit: */ 0
                                 , /* upper_limit: */ problem.getINF()
                                 , /*      coeffs: */ &cluster_constraint );
            for ( LidT j = 0; j != quadConstraints.size(); ++j )
                problem.addQConstraint( constraintId, quadConstraints[j].row(), quadConstraints[j].col(), quadConstraints[j].value() );
        }

    } //...dId constraint
    std::cout << "[" << __func__ << "]: " << "lvl2 constr end..." << std::endl; fflush(stdout);

    std::cout << "[" << __func__ << "]: " << "lvl2 pw start..." << std::endl; fflush(stdout);
    // ____________________________________________________
    // dId pw cost
    {
        //#pragma omp parallel for num_threads(GF2_MAX_OMP_THREADS)
        for ( auto it0 = dIdsVarIds.begin(); it0 != dIdsVarIds.end(); ++it0 )
            for ( auto it1 = dIdsVarIds.begin(); it1 != dIdsVarIds.end(); ++it1 )
            {
                // first: dId
                // second: varId

                if ( it1->second == it0->second ) continue; // self pw cost is 0

                _PrimitiveT const* p0 = dIdsPrims[it0->first];
                _PrimitiveT const* p1 = dIdsPrims[it1->first];
                //
                //_Scalar score = weights(1) * sqrt( GF2::angleInRad(p0->template dir(), p1->template dir()) );
                //_Scalar score = weights(1) * sqrt( MyPrimitivePrimitiveAngleFunctor::eval( *p0, *p1, angles ) ); //changed 16:11 15/01/2015
                _Scalar score = weights(1) * calcPwCost<_Scalar>( *p0, *p1, angles );
//#               pragma omp critical (PS_PROBLEM)
                {
                    problem.addQObjective( it0->second, it1->second
                                           , score // score
                                           );
                }
            }
    } //...dId pw cost
    std::cout << "[" << __func__ << "]: " << "lvl2 pw end..." << std::endl; fflush(stdout);

    std::cout << "[" << __func__ << "]: " << "init solution..." << std::endl; fflush(stdout);
    // ____________________________________________________
    // Initial solution
    {
        if ( chosen_varids.size() )
        {
            OptProblemT::SparseMatrix x0( problem.getVarCount(), 1 );
            for ( std::set<LidT>::const_iterator it = chosen_varids.begin(); it != chosen_varids.end(); ++it )
            {
                x0.insert( *it, 0 ) = 1;
            }
            problem.setStartingPoint( x0 );
        }
    } //...Initial solution
    std::cout << "[" << __func__ << "]: " << "init solution done" << std::endl; fflush(stdout);

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

    //! \brief              Adds constraints to \p problem so, that each patch (prims[i] that have the same _PrimitiveT::TAGS::GID) has at least one member j (prims[i][j]) selected.
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

                if ( verbose && (lid1 == 0) ) std::cout << "[" << __func__ << "]: " << "Constraining " << prims[lid][lid1].getTag( _PrimitiveT::TAGS::GID ) << " to choose one of ";
                coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = 1.0;
                non_zero_line = true;
                if ( verbose ) std::cout << prims[lid][lid1].getTag( _PrimitiveT::TAGS::DIR_GID ) << ", ";
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

    //____________________________________________DataCosts____________________________________________

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
                            , _Scalar              const  scale
                            , _Scalar              const freq_weight
                            , bool                 const verbose )
    {
        typedef typename _AssocT::key_type IntPair;

        const _Scalar w_mod_base     = ProblemSetupParams<_Scalar>::w_mod_base,
                      w_mod_base_inv = _Scalar(1.) - w_mod_base;

        //int err = EXIT_SUCCESS;

        // count patch populations
        GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        processing::getPopulations( populations, points );

        // INSTANCES
        // added 17/09/2014 by Aron
        // modified 22/09/2014 by Aron
#if 0
        std::map< GidT, LidT > dir_instances;
        LidT active_count = 0;
        for ( size_t lid = 0; lid != prims.size(); ++lid )
            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
            {
                if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) != _PrimitiveT::STATUS_VALUES::ACTIVE )
                    continue;

                ++dir_instances[ prims[lid][lid1].getTag(_PrimitiveT::TAGS::DIR_GID) ];
                ++active_count;
            }
#endif

#       pragma omp parallel for num_threads(GF2_MAX_OMP_THREADS)
        for ( size_t lid = 0; lid < prims.size(); ++lid )
        {
            // check, if any directions for patch
            if ( !prims[lid].size() )
            {
                //std::cerr << "[" << __func__ << "]: " << "no directions for patch[" << lid << "]! This shouldn't happen. Skipping patch..." << std::endl;
                continue;
            }

            // cache patch group id to match with point group ids
            const GidT gid = prims[lid][0].getTag( _PrimitiveT::TAGS::GID );

            // for each direction
            for ( size_t lid1 = 0; lid1 < prims[lid].size(); ++lid1 )
            {
                if ( prims[lid][lid1].getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                    continue;

                typename _PrimitiveT::ExtremaT extrema;
                int err = prims[lid][lid1].template getExtent<_PointPrimitiveT>
                        ( extrema
                        , points
                        , scale
                        , populations[gid].size() ? &(populations[gid]) : NULL );

                // point count for normalization
                unsigned cnt = 0;
                // data-cost coefficient (output)
                _Scalar unary_i = _Scalar(0);
                // for each point, check if assigned to main patch (TODO: move assignment test to earlier, it's indep of lid1)
                for ( size_t pid = 0; pid != points.size(); ++pid )
                {
                    //if ( points[pid].getTag( PointPrimitiveT::TAGS::GID ) == static_cast<int>(lid) )
                    if ( points[pid].getTag( _PointPrimitiveT::TAGS::GID ) == gid )
                    {
                        // if within scale, add unary cost

                        // changed by Aron on 6/1/2015
                         _Scalar dist = std::numeric_limits<_Scalar>::max();
                        if ( err == EXIT_SUCCESS )
                        {
                            dist = MyPointFiniteLineDistanceFunctor::eval( extrema, prims[lid][lid1], points[pid].template pos() );
                        }
                        else
                        {
                            //dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );

                            dist = 2.; // we don't want an empty primitive
                        }
                        unary_i += dist * dist; //changed on 18/09/14
                        ++cnt;              // normalizer
                    }
                } // for points

                // average data cost
                _Scalar coeff = cnt ? /* unary: */ weights(0) * unary_i / _Scalar(cnt)
                                    : /* unary: */ weights(0) * _Scalar(2);            // add large weight, if no points assigned
//                std::cout << "coeff(" << gid
//                          << ","
//                          << prims[lid][lid1].getTag(_PrimitiveT::TAGS::DIR_GID) << ") = "
//                          << coeff
//                          << ", unary: " << unary_i << ", cnt: " << cnt << ", weights(0): " << weights(0)
//                          << std::endl;

#if 0
                // prefer dominant directions
                if ( freq_weight > _Scalar(0.) )
                {
                    const DidT dir_gid = prims[lid][lid1].getTag( _PrimitiveT::TAGS::GID );

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
#endif

                // complexity cost:
                coeff += weights(2); // changed by Aron on 21/9/2014

                // add to problem
#               pragma omp critical (PS_PROBLEM)
                {
                    problem.addLinObjective( /* var_id: */ lids_varids.at( IntPair(lid,lid1) )
                                           , /*  value: */ coeff );
                }
            } //...for each direction
        } //...for each patch

        return EXIT_SUCCESS;
    } //...associationBasedDataCost


} //...namespace ProblemSetup
} //...namespace GF2

#endif // GF2_PROBLEMSETUP_HPP
