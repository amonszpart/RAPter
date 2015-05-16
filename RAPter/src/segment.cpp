#include "rapter/typedefs.h"                    // _2d::, _3d::, Scalar, PointPrimitiveT, PointContainerT
#include "rapter/optimization/segmentation.h"   // segmentCli, orientPoints, patchify
#include "rapter/util/parse.h"                  // console::

int segment( int argc, char**argv )
{
    if ( rapter::console::find_switch(argc,argv,"--segment3D") )
    {
        return rapter::Segmentation::segmentCli< rapter::_3d::PrimitiveT
                                            , rapter::_3d::PrimitiveContainerT
                                            , rapter::PointPrimitiveT
                                            , rapter::PointContainerT
                                            , rapter::Scalar>
                                            ( argc, argv );
    }
    else
    {
        return rapter::Segmentation::segmentCli< rapter::_2d::PrimitiveT
                                            , rapter::_2d::PrimitiveContainerT
                                            , rapter::PointPrimitiveT
                                            , rapter::PointContainerT
                                            , rapter::Scalar>
                                            ( argc, argv );
    } //...find_switch
} //...segment()
