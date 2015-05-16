#include "globfit2/globOpt_types.h"
//#include "rapter/primitives/impl/planePrimitive.hpp"
#include "globopt/evaluation/impl/assignPointsToTriangles.hpp"

namespace rapter
{
    template
    int assignPointsToTriangles< rapter::_3d::PrimitiveVectorT
                                 , rapter::_3d::PrimitiveMapT
                                 , rapter::PointContainerT
                                 , rapter::PclCloudT>( int argc, char** argv );
} //...ns rapter
