#include "globfit2/optimization/segmentation.h" // segmentCli, orientPoints, patchify

#include "globfit2/globOpt_types.h" // _2d::, _3d::, Scalar, PointPrimitiveT, PointContainerT
#include "globfit2/util/parse.h" // console::

int segment( int argc, char**argv )
{
    if ( GF2::console::find_switch(argc,argv,"--segment3D") )
    {
//    typedef GF2::PlanePrimitive              PrimitiveT;
//    typedef std::vector< PrimitiveT >        InnerContainerT;
//    typedef std::vector< InnerContainerT >   PrimitiveContainerT; // TODO: map
//    typedef GF2::PointPrimitive              PointPrimitiveT;
//    typedef std::vector< PointPrimitiveT >   PointContainerT;
//    typedef PrimitiveT::Scalar               Scalar;

        return GF2::Segmentation::segmentCli< GF2::_3d::PrimitiveT
                                            , GF2::_3d::PrimitiveContainerT
                                            , GF2::PointPrimitiveT
                                            , GF2::PointContainerT
                                            , GF2::Scalar>
                                            ( argc, argv );
    }
    else
    {

//    typedef GF2::LinePrimitive2              PrimitiveT;
//    typedef std::vector< PrimitiveT >        InnerContainerT;
//    typedef std::vector< InnerContainerT > PrimitiveContainerT; // TODO: map
//    typedef GF2::PointPrimitive              PointPrimitiveT;
//    typedef std::vector< PointPrimitiveT >   PointContainerT;
//    typedef PrimitiveT::Scalar               Scalar;
        return GF2::Segmentation::segmentCli< GF2::_2d::PrimitiveT
                                            , GF2::_2d::PrimitiveContainerT
                                            , GF2::PointPrimitiveT
                                            , GF2::PointContainerT
                                            , GF2::Scalar>
                                            ( argc, argv );
    } //...find_switch
} //...segment()
