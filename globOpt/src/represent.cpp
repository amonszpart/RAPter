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

    GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
    processing::getPopulations( populations, points );

    std::cout << "gids: " << patches.size() << ", points: " << points.size() << ", scale: " << params.scale << std::endl;

    SpatialSqrtPrimitivePrimitiveEnergyFunctor<_FiniteFiniteDistFunctor, _PointContainerT, Scalar,_PrimitiveT>
            ppDistFunctor( params.angles, points, params.scale );

    graph::EdgeListT<Scalar> edgesList;

    for ( typename PrimitiveMapT::Iterator it0(patches); it0.hasNext(); it0.step() )
        for ( typename PrimitiveMapT::Iterator it1(patches); it1.hasNext(); it1.step() )
    {
        if ( it1 < it0 ) continue; // WARNING: hack, assumes gids are coming sorted increasingly

        const DidT did0 = it0->getTag(_PrimitiveT::DIR_GID);
        const DidT did1 = it1->getTag(_PrimitiveT::DIR_GID);

        std::cout << "<" << it0.getGid() << "," << it0.getLid() << "," << did0 << "," << it0.getUniqueId() << ">"
                  << "<" << it1.getGid() << "," << it1.getLid() << "," << did1 << "," << it1.getUniqueId() << ">"
                  << std::endl;


        ExtentsT extrema, extrema2;
        int err = it0->template getExtent<_PointPrimitiveT>( extrema
                                                          , points
                                                          , params.scale
                                                          , &(populations[it0.getGid()]) )
                + it1->template getExtent<_PointPrimitiveT>( extrema2
                                                            , points
                                                            , params.scale
                                                            , &(populations[it1.getGid()]) );
        if ( EXIT_SUCCESS == err )
        {
            Scalar invDist = ppDistFunctor.evalSpatial( *it0, extrema
                                                      , *it1, extrema2 );
            if ( invDist > Scalar(0.) && (did0 == did1) )
            {
                edgesList.insert( EdgeT(it0.getUniqueId(), it1.getUniqueId(), invDist) );
            }
        }
    }

    GraphT graph( edgesList );
    graph.draw( "representGraph.gv", /* show: */ true );

    typename GraphT::ClustersT clusters;
    graph.getClusters( clusters, 2 );
    graph.showClusters( clusters, "representClusters.gv", /* show: */ true );



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


