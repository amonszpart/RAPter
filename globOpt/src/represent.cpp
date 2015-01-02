#include "globfit2/io/inputParser.hpp"
#include "globfit2/util/parse.h"
#include "globfit2/globOpt_types.h"
#include "globfit2/my_types.h"
#include "globfit2/processing/graph.hpp"
#include "globfit2/util/containers.hpp" // class PrimitiveContainer

namespace GF2
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
    PointContainerT         points;
    PclCloudPtrT            pcl_cloud;
    _PrimitiveContainerT    prims;
    PrimitiveMapT           patches;
    RepresentParams<Scalar> params;
    typedef typename _PrimitiveT::ExtentsT                   ExtentsT;

    // Graphs
    typedef Graph<Scalar, typename MyGraphConfig<Scalar>::UndirectedGraph > GraphT;
    typedef typename graph::EdgeT<Scalar>                                   EdgeT;
    typedef typename GraphT::ComponentListT                                 ComponentListT;

    bool valid_input = GF2::parseInput<InnerPrimitiveContainerT,PclCloudT>(
                        points, pcl_cloud, prims, patches, params, argc, argv );
    AnglesT angle_gens;
    parseAngles( params.angles, argc, argv, &angle_gens );

    GidPidVectorMap populations;
    processing::getPopulations( populations, points );

    std::cout << "[" << __func__ << "]: " << "gids: " << patches.size() << ", points: " << points.size() << ", scale: " << params.scale << std::endl;

    // finite distance calculation between patches
    SpatialSqrtPrimitivePrimitiveEnergyFunctor<_FiniteFiniteDistFunctor, _PointContainerT, Scalar,_PrimitiveT>
            ppDistFunctor( params.angles, points, params.scale );

    // create spatial+directionId clustered graph via this edge list
    graph::EdgeListT<Scalar> edgesList;

    // enumerate all unique primitive pairs
    for ( typename PrimitiveMapT::Iterator it0(patches); it0.hasNext(); it0.step() )
        for ( typename PrimitiveMapT::Iterator it1(patches); it1.hasNext(); it1.step() )
    {
        // skip reverse pairs (0-1 allows 1-0 to be skipped)
        if ( it1 < it0 ) continue; // WARNING: hack, assumes gids are coming sorted increasingly

        // cache direction ids
        const DidT did0 = it0->getTag(_PrimitiveT::DIR_GID);
        const DidT did1 = it1->getTag(_PrimitiveT::DIR_GID);

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

        // add edge, if same direction id (colour) and close to eachother
        if ( EXIT_SUCCESS == err )
        {
            // Originally, evalSpatial returns 1 if they are at the same spot, and 0 if they are 2 x scale away.
            Scalar invDist = ppDistFunctor.evalSpatial( *it0, extrema
                                                      , *it1, extrema1 );
            if ( (did0 == did1) && invDist > Scalar(0.) ) // same colour and closer than 2xscale
                edgesList.insert( EdgeT(it0.getUniqueId(), it1.getUniqueId(), /* not used: */ invDist) );
        } // if extrema exist
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

    // WORK


    return !valid_input;
} //...representCli

} //...GF2

int represent( int argc, char** argv )
{
#if 1
     if ( GF2::console::find_switch(argc,argv,"--represent3D") )
     {
         return GF2::representCli< GF2::_3d::PrimitiveContainerT
                 , GF2::PointContainerT
                 , GF2::_3d::PrimitiveT
                 , GF2::PointPrimitiveT
                 , GF2::_3d::MyFinitePlaneToFinitePlaneCompatFunctor
                 >( argc, argv );
     }
     else
     {
         return GF2::representCli< GF2::_2d::PrimitiveContainerT
                                 , GF2::PointContainerT
                                 , GF2::_2d::PrimitiveT
                                 , GF2::PointPrimitiveT
                                 , GF2::_2d::MyFiniteLineToFiniteLineCompatFunctor
                                 >( argc, argv );
     }
#endif

     return 0;
}


