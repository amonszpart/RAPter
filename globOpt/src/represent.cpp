#include "rapter/io/inputParser.hpp"
#include "rapter/util/parse.h"
#include "rapter/globOpt_types.h"
#include "rapter/my_types.h"
#include "rapter/processing/graph.hpp"
#include "rapter/util/containers.hpp" // class PrimitiveContainer

namespace rapter
{

template <typename _Scalar>
struct RepresentParams
{
    _Scalar scale;
    AnglesT angles;
};

template <
           class _PrimitiveContainerT
         , class _PointContainerT
         , class _PrimitiveT /*         = typename _PrimitiveContainerT::value_type::value_type*/
         , class _PointPrimitiveT     /*= typename _PointContainerT::value_type*/
         , class _FiniteFiniteDistFunctor
         >
static inline int representCli( int argc, char** argv )
{
    // input
    typedef typename _PrimitiveContainerT::value_type        InnerPrimitiveContainerT;
    typedef containers::PrimitiveContainer<_PrimitiveT>      PrimitiveMapT;
    typedef typename _PrimitiveT::Scalar                     Scalar;
    typedef typename _PrimitiveT::ExtremaT                   ExtremaT;
    PointContainerT         points;
    PclCloudPtrT            pcl_cloud;
    _PrimitiveContainerT    prims;
    PrimitiveMapT           patches;
    RepresentParams<Scalar> params;

    // Graphs
    typedef Graph<Scalar, typename MyGraphConfig<Scalar>::UndirectedGraph > GraphT;
    typedef typename graph::EdgeT<Scalar>                                   EdgeT;
    typedef typename GraphT::ComponentListT                                 ComponentListT;
    typedef typename GraphT::ClustersT                                      ClustersT;

    int ret = rapter::parseInput<InnerPrimitiveContainerT,PclCloudT>(
                points, pcl_cloud, prims, patches, params, argc, argv );
    std::cout << "[" << __func__ << "]: " << "parseInput ret: " << ret << std::endl;
    bool valid_input = (EXIT_SUCCESS == ret);

    AnglesT angle_gens;
    valid_input &= (EXIT_SUCCESS == parseAngles(params.angles, argc, argv, &angle_gens) );

    if ( !valid_input )
    { std::cout << "Usage: --represent[3D] -p prims.csv -a points_primitives.csv -sc scale --cloud cloud.ply --angle-gens 90" << std::endl; return EXIT_FAILURE; }

    GidPidVectorMap populations;
    processing::getPopulations( populations, points );

    std::cout << "[" << __func__ << "]: " << "gids: " << patches.size() << ", points: " << points.size() << ", scale: " << params.scale << std::endl;

#if 0
    // finite distance calculation between patches
    SpatialSqrtPrimitivePrimitiveEnergyFunctor<_FiniteFiniteDistFunctor, _PointContainerT, Scalar,_PrimitiveT>
            ppDistFunctor( params.angles, points, params.scale );

    // create spatial+directionId clustered graph via this edge list
    graph::EdgeListT<Scalar> edgesList;
    std::map< UidT, _PrimitiveT const* > primitivesByUid;

    // enumerate all unique primitive pairs
    for ( typename PrimitiveMapT::Iterator it0(patches); it0.hasNext(); it0.step() )
        for ( typename PrimitiveMapT::Iterator it1(patches); it1.hasNext(); it1.step() )
    {
        // skip reverse pairs (0-1 allows 1-0 to be skipped)
        if ( it1 < it0 ) continue; // WARNING: hack, assumes gids are coming sorted increasingly

        // cache direction ids
        const DidT did0 = it0->getTag(_PrimitiveT::TAGS::DIR_GID);
        const DidT did1 = it1->getTag(_PrimitiveT::TAGS::DIR_GID);

        // log
        std::cout << "<" << it0.getGid() << "," << it0.getLid() << "," << did0 << "," << it0.getUniqueId() << ">"
                  << "<" << it1.getGid() << "," << it1.getLid() << "," << did1 << "," << it1.getUniqueId() << ">"
                  << std::endl;

        // get spatial extent (is usually cached in the function)
        ExtentsT extrema, extrema1;
        int err = it0->template getExtent<_PointPrimitiveT>( extrema
                                                           , points
                                                           , params.scale
                                                           , &(populations[it0.getGid()]) )
                + it1->template getExtent<_PointPrimitiveT>( extrema1
                                                           , points
                                                           , params.scale
                                                           , &(populations[it1.getGid()]) );

        // add edge, if same direction id (colour) and close to each other
        if ( EXIT_SUCCESS == err )
        {
            // Originally, evalSpatial returns 1 if they are at the same spot, and 0 if they are 2 x scale away.
            Scalar invDist = ppDistFunctor.evalSpatial( *it0, extrema
                                                      , *it1, extrema1 );
            if ( (did0 == did1) && invDist > Scalar(0.) ) // same colour and closer than 2xscale
                edgesList.insert( EdgeT(it0.getUniqueId(), it1.getUniqueId(), /* not used: */ invDist) );

            if ( primitivesByUid.find(it0.getUniqueId()) == primitivesByUid.end() )
                primitivesByUid[ it0.getUniqueId() ] = &(*it0);
            if ( primitivesByUid.find(it1.getUniqueId()) == primitivesByUid.end() )
                primitivesByUid[ it1.getUniqueId() ] = &(*it1);
        } //...if extrema exist
    } //...all primitive pairs

    // build fixed size graph from edge list
    GraphT graph( edgesList );
    // plot (debug)
    graph.draw( "representGraph.gv", /* show: */ true );

    // extract connected components with at least 2 nodes
    typename GraphT::ClustersT clusters;
    graph.getClusters( clusters, /* minimum primitive count: */ 2 );
    // plot (debug)
    graph.showClusters( clusters, "representClusters.gv", /* show: */ true );


    // for each cluster, select the largest (for now) representative
    for ( typename ClustersT::const_iterator clustersIt = clusters.begin();
          clustersIt != clusters.end(); ++clustersIt )
    {
        // first: clusterId
        // second: vector<nodeId>

        typedef typename ClustersT::mapped_type InnerT;
        InnerT const& nodes = clustersIt->second;
        for ( typename InnerT::const_iterator nodeIt = nodes.begin();
              nodeIt != nodes.end(); ++nodeIt )
        {
            // nodeIt is a uniqueId in the graph
            std::cout << "fetching " << *nodeIt;
            _PrimitiveT const* prim = primitivesByUid.at( *nodeIt );
            std::cout << "OK" << std::endl;
            fflush(stdout);
        }
    }
#endif
    /// WORK
    // output representatives
    _PrimitiveContainerT outPrims;
    std::set<GidT> activeGids;

    typedef typename Eigen::Matrix<Scalar,Eigen::Dynamic,1> SpatialSignifT;
    SpatialSignifT spatialSignif(1,1); // tmp
    typedef typename std::pair<SpatialSignifT,_PrimitiveT const*> SizedPrimT;
    std::map< DidT, SizedPrimT > maxSpatialSignifs; // "size"

    for ( typename PrimitiveMapT::Iterator it0(patches); it0.hasNext(); it0.step() )
    {
        if ( it0->getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL ) continue; // added 9 / 1 / 2015

        // cache
        _PrimitiveT const* prim = &(*it0);
        const DidT did = prim->getTag( _PrimitiveT::TAGS::DIR_GID );
        // calc size
        prim->getSpatialSignificance( spatialSignif, points, params.scale, &(populations[prim->getTag(_PrimitiveT::TAGS::GID)]) );

        // insert, if did unseen
        if ( maxSpatialSignifs.find( did ) == maxSpatialSignifs.end() )
            maxSpatialSignifs[ did ] = SizedPrimT( spatialSignif, prim );
        // or replace max, if larger
        else if ( spatialSignif(0) > maxSpatialSignifs[did].first(0) )
        {
            maxSpatialSignifs[did].first  = spatialSignif; // store size
            maxSpatialSignifs[did].second = prim;          // store primitive
        } //...if bigger
    } //...all primitives

    // record to output
    for ( auto it = maxSpatialSignifs.begin(); it != maxSpatialSignifs.end(); ++it )
    {
        //  first: did
        // second: SizedPrimT (<size,prim>)
        _PrimitiveT const* prim = it->second.second;
        containers::add( outPrims, prim->getTag(_PrimitiveT::TAGS::GID), *(prim) );
        activeGids.insert( prim->getTag(_PrimitiveT::TAGS::GID) );
    }

    std::cout << "[" << __func__ << "]: " << "saved to representatives.csv" << std::endl;
    io::savePrimitives<_PrimitiveT,typename InnerPrimitiveContainerT::const_iterator>( outPrims, "representatives.csv" );

    // ___POINTS___
    _PointContainerT outPoints( points ); // need reassignment
    // clear all points' assignment that don't have selected primitives
    for ( auto it = outPoints.begin(); it != outPoints.end(); ++it )
    {
        if ( activeGids.find( it->getTag(_PointPrimitiveT::TAGS::GID) ) == activeGids.end() )
            it->setTag( _PointPrimitiveT::TAGS::GID, _PointPrimitiveT::LONG_VALUES::UNSET );
    }
    std::string assocPath( "points_representatives.csv" );
    io::writeAssociations<_PointPrimitiveT>( outPoints, assocPath );

    return !valid_input;
} //...representCli

/*! \brief Takes ran representatives, and changes prims accordingly.
 */
template <
           class _PrimitiveContainerT
         , class _PointContainerT
         , class _PrimitiveT /*         = typename _PrimitiveContainerT::value_type::value_type*/
         , class _PointPrimitiveT     /*= typename _PointContainerT::value_type*/
         , class _FiniteFiniteDistFunctor
         >
static inline int representBackCli( int argc, char** argv )
{
    const int verbose = 1;
    // input
    typedef typename _PrimitiveContainerT::value_type        InnerPrimitiveContainerT;
    typedef containers::PrimitiveContainer<_PrimitiveT>      PrimitiveMapT;
    typedef typename _PrimitiveT::Scalar                     Scalar;
    typedef typename _PrimitiveT::ExtremaT                   ExtremaT;
    PointContainerT         points;
    PclCloudPtrT            pcl_cloud;
    PrimitiveMapT           patches,reprPatches;
    RepresentParams<Scalar> params;

    // Graphs
    typedef Graph<Scalar, typename MyGraphConfig<Scalar>::UndirectedGraph > GraphT;
    typedef typename graph::EdgeT<Scalar>                                   EdgeT;
    typedef typename GraphT::ComponentListT                                 ComponentListT;
    typedef typename GraphT::ClustersT                                      ClustersT;

    // read input
    bool valid_input = true;
    {
        _PrimitiveContainerT prims;
        valid_input = ( EXIT_SUCCESS == rapter::parseInput<InnerPrimitiveContainerT,PclCloudT>(
                                 points, pcl_cloud, prims, patches, params, argc, argv, true ) );
        if ( !valid_input ) std::cout << "failed first" << std::endl;
    }

    // read angles
    AnglesT angle_gens;
    valid_input &= (EXIT_SUCCESS == parseAngles(params.angles, argc, argv, &angle_gens) );
    if ( !valid_input ) std::cout << "failed angles" << std::endl;

    // read finished representatives
    std::string representativesPath;
    valid_input &= (0 <= rapter::console::parse_argument(argc, argv, "--repr", representativesPath));
    if ( !valid_input ) std::cout << "failed reprpath" << std::endl;
    {
        _PrimitiveContainerT representatives;
        valid_input &= (EXIT_SUCCESS == rapter::io::readPrimitives<_PrimitiveT, InnerPrimitiveContainerT>( representatives, representativesPath, &reprPatches) );
        if ( !valid_input ) std::cout << "failed read" << std::endl;
    }

    if ( !valid_input )
    { std::cerr << "Usage: --representBack[3D] -p prims.csv -a points_primitives.csv -sc scale --cloud cloud.ply --repr representatives.bonmin.csv --angle-gens 90" << std::endl; return EXIT_FAILURE; }

    GidPidVectorMap populations;
    processing::getPopulations( populations, points );

    std::cout << "reprback" << std::endl;

    // loop over representative output to check for ID matches
    std::map<DidT, _PrimitiveT const*> subs;
    for ( typename PrimitiveMapT::Iterator it(reprPatches); it.hasNext(); it.step() )
    {
        if ( patches[ it.getGid() ].size() > 1 )
        {
            std::cerr << "[" << __func__ << "]: " << "can't handle more than one primitive per patch( gid: " << it.getGid() << "!" << std::endl;
            //return 1;
        }
        //for ( int )
        const int dId = patches[ it.getGid() ].at(0).getTag(_PrimitiveT::TAGS::DIR_GID);
        std::cout << "did at " << it.getGid() << " is " << dId << std::endl;

        if ( it.getDid() != dId )
        {
            if ( verbose ) std::cout << "found subst: " << it.getDid() << " instead of " << dId << " at gid " << it.getGid() << std::endl;
            if (    (subs.find( dId )                              != subs.end() )
                 && (subs[dId]->getTag(_PrimitiveT::TAGS::DIR_GID) != it.getDid()) )
                std::cerr << "duplicate subs for patch " << it.getGid() << ": " << subs[dId]->getTag(_PrimitiveT::TAGS::DIR_GID) << ", " << it.getDid() << "!" << std::endl;
            else
                subs[ dId ] = &(*it);
        }
    }

    _PrimitiveContainerT outPrims;
    for ( typename PrimitiveMapT::Iterator it(patches); it.hasNext(); it.step() )
    {
        bool copy = false;

        copy = (it->getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL );

        if ( !copy && (subs.find(it.getDid()) != subs.end()) )
        {
            _PrimitiveT const* subExample = subs[ it.getDid() ]; // pattern, to copy direction from
            std::cout << "substituting " << it.getGid() << "," << it.getDid() << " <- " << subExample->getTag(_PrimitiveT::TAGS::DIR_GID) << std::endl;

            int closest_angle_id = 0;
            _PrimitiveT sub;
            Scalar angdiff = MyPrimitivePrimitiveAngleFunctor::template eval<Scalar>( *it, *subExample
                                                                                    , params.angles
                                                                                    , &closest_angle_id );
            if ( angdiff > 0.1 )
                std::cerr << "[" << __func__ << "]: " << "really? substitute dId, but best angdiff: " << angdiff << std::endl;

            if ( !it->generateFrom(/* out */ sub, /* example: */ *subExample, closest_angle_id, params.angles, /* sign, unused: */ Scalar(1.)) )
            {
                std::cerr << "[" << __func__ << "]: " << "could not generate candidate" << std::endl;
                copy = true;
            }
            else
            {
                sub.setTag( _PrimitiveT::TAGS::STATUS, it->getTag(_PrimitiveT::TAGS::STATUS) );
                containers::add( outPrims, it.getGid(), sub );
            }
        }
        else
            copy = true;

        if ( copy )
        {
            containers::add( outPrims, it.getGid(), *it );
        }
    }

    return io::savePrimitives<_PrimitiveT, typename InnerPrimitiveContainerT::const_iterator >( outPrims, "subs.csv" );
} //...representBack

} //...ns rapter

int represent( int argc, char** argv )
{
#if 1
     if ( rapter::console::find_switch(argc,argv,"--represent3D") )
     {
         return rapter::representCli< rapter::_3d::PrimitiveContainerT
                 , rapter::PointContainerT
                 , rapter::_3d::PrimitiveT
                 , rapter::PointPrimitiveT
                 , rapter::_3d::MyFinitePlaneToFinitePlaneCompatFunctor
                 >( argc, argv );
     }
     else if ( rapter::console::find_switch(argc,argv,"--represent") )
     {
         return rapter::representCli< rapter::_2d::PrimitiveContainerT
                                 , rapter::PointContainerT
                                 , rapter::_2d::PrimitiveT
                                 , rapter::PointPrimitiveT
                                 , rapter::_2d::MyFiniteLineToFiniteLineCompatFunctor
                                 >( argc, argv );
     }
     else if ( rapter::console::find_switch(argc,argv,"--representBack3D") )
     {
         return rapter::representBackCli< rapter::_3d::PrimitiveContainerT
                 , rapter::PointContainerT
                 , rapter::_3d::PrimitiveT
                 , rapter::PointPrimitiveT
                 , rapter::_3d::MyFinitePlaneToFinitePlaneCompatFunctor
                 >( argc, argv );
     }
     else if ( rapter::console::find_switch(argc,argv,"--representBack") )
     {
         return rapter::representBackCli< rapter::_2d::PrimitiveContainerT
                 , rapter::PointContainerT
                 , rapter::_2d::PrimitiveT
                 , rapter::PointPrimitiveT
                 , rapter::_2d::MyFiniteLineToFiniteLineCompatFunctor
                 >( argc, argv );
     }
     else
         std::cerr << "[" << __func__ << "]: " << "switch error" << std::endl;

#endif

     return 0;
}


