#include "rapter/optimization/segmentation.h" // segmentCli, orientPoints, patchify

#include "rapter/globOpt_types.h" // _2d::, _3d::, Scalar, PointPrimitiveT, PointContainerT
#include "rapter/util/parse.h" // console::

int segment( int argc, char**argv )
{
    if ( rapter::console::find_switch(argc,argv,"--segment3D") )
    {
//    typedef rapter::PlanePrimitive              PrimitiveT;
//    typedef std::vector< PrimitiveT >        InnerContainerT;
//    typedef std::vector< InnerContainerT >   PrimitiveContainerT; // TODO: map
//    typedef rapter::PointPrimitive              PointPrimitiveT;
//    typedef std::vector< PointPrimitiveT >   PointContainerT;
//    typedef PrimitiveT::Scalar               Scalar;

        return rapter::Segmentation::segmentCli< rapter::_3d::PrimitiveT
                                            , rapter::_3d::PrimitiveContainerT
                                            , rapter::PointPrimitiveT
                                            , rapter::PointContainerT
                                            , rapter::Scalar>
                                            ( argc, argv );
    }
    else
    {

//    typedef rapter::LinePrimitive2              PrimitiveT;
//    typedef std::vector< PrimitiveT >        InnerContainerT;
//    typedef std::vector< InnerContainerT > PrimitiveContainerT; // TODO: map
//    typedef rapter::PointPrimitive              PointPrimitiveT;
//    typedef std::vector< PointPrimitiveT >   PointContainerT;
//    typedef PrimitiveT::Scalar               Scalar;
        return rapter::Segmentation::segmentCli< rapter::_2d::PrimitiveT
                                            , rapter::_2d::PrimitiveContainerT
                                            , rapter::PointPrimitiveT
                                            , rapter::PointContainerT
                                            , rapter::Scalar>
                                            ( argc, argv );
    } //...find_switch
} //...segment()
