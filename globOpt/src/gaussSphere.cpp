#include "globfit2/io/inputParser.hpp"
#include "globfit2/util/parse.h"
#include "globfit2/globOpt_types.h"
#include "globfit2/my_types.h"
#include "globfit2/processing/graph.hpp"
#include "globfit2/util/containers.hpp" // class PrimitiveContainer

namespace GF2
{

template <typename _Scalar>
struct GaussSphereParams
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
static inline int gaussSphereCli( int argc, char** argv )
{
    // input
    typedef typename _PrimitiveContainerT::value_type        InnerPrimitiveContainerT;
    typedef containers::PrimitiveContainer<_PrimitiveT>      PrimitiveMapT;
    typedef typename _PrimitiveT::Scalar                     Scalar;
    typedef typename _PrimitiveT::ExtentsT                   ExtentsT;
    PointContainerT         points;
    PclCloudPtrT            pcl_cloud;
    _PrimitiveContainerT    prims;
    PrimitiveMapT           patches;
    GaussSphereParams<Scalar> params;

    // Graphs
    typedef Graph<Scalar, typename MyGraphConfig<Scalar>::UndirectedGraph > GraphT;
    typedef typename graph::EdgeT<Scalar>                                   EdgeT;
    typedef typename GraphT::ComponentListT                                 ComponentListT;
    typedef typename GraphT::ClustersT                                      ClustersT;

    int ret = GF2::parseInput<InnerPrimitiveContainerT,PclCloudT>(
                points, pcl_cloud, prims, patches, params, argc, argv );
    std::cout << "[" << __func__ << "]: " << "parseInput ret: " << ret << std::endl;
    bool valid_input = (EXIT_SUCCESS == ret);

    AnglesT angle_gens;
    valid_input &= (EXIT_SUCCESS == parseAngles(params.angles, argc, argv, &angle_gens) );

    if ( !valid_input )
    { std::cout << "Usage: [--3D] -p prims.csv -a points_primitives.csv -sc scale --cloud cloud.ply --angle-gens 90" << std::endl; return EXIT_FAILURE; }

    GidPidVectorMap populations;
    processing::getPopulations( populations, points );

    std::cout << "[" << __func__ << "]: " << "gids: " << patches.size() << ", points: " << points.size() << ", scale: " << params.scale << std::endl;

    //for ( )

    return EXIT_SUCCESS;
} //...gaussSphereCli
} //...ns GF2

int main( int argc, char *argv[] )
{
    if ( GF2::console::find_switch(argc,argv,"--3D") || GF2::console::find_switch(argc,argv,"--3d") )
    {
        return GF2::gaussSphereCli< GF2::_3d::PrimitiveContainerT
                , GF2::PointContainerT
                , GF2::_3d::PrimitiveT
                , GF2::PointPrimitiveT
                , GF2::_3d::MyFinitePlaneToFinitePlaneCompatFunctor
                >( argc, argv );
    }
    else
    {
        return GF2::gaussSphereCli< GF2::_2d::PrimitiveContainerT
                                , GF2::PointContainerT
                                , GF2::_2d::PrimitiveT
                                , GF2::PointPrimitiveT
                                , GF2::_2d::MyFiniteLineToFiniteLineCompatFunctor
                                >( argc, argv );
    }

    return 0;
}

