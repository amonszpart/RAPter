#include "rapter/optimization/merging.h"  // mergeCli

#include "rapter/globOpt_types.h"         // _2d::PrimitiveT, _3d::PrimitiveT
#include "rapter/util/parse.h"            // find_switch
#include "rapter/primitives/impl/planePrimitive.hpp"

int merge( int argc, char** argv )
{
    if ( rapter::console::find_switch(argc,argv,"--merge3D") )
    {
        return rapter::Merging::mergeCli< rapter::_3d::PrimitiveContainerT
                                     , rapter::PointContainerT
                                     , rapter::Scalar
                                     , rapter::PointPrimitiveT
                                     , rapter::_3d::PrimitiveT
                                     >( argc, argv );
    }
    else
    {
        return rapter::Merging::mergeCli< rapter::_2d::PrimitiveContainerT
                                     , rapter::PointContainerT
                                     , rapter::Scalar
                                     , rapter::PointPrimitiveT
                                     , rapter::_2d::PrimitiveT
                                     >( argc, argv );
    } //...if find_switch
} //...merge()

