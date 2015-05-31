#include "rapter/typedefs.h"
//#include "rapter/primitives/impl/planePrimitive.hpp"
#include "rapter/evaluation/impl/assignPointsToTriangles.hpp"

namespace rapter
{
    template
    int assignPointsToTriangles< rapter::_3d::PrimitiveVectorT
                                 , rapter::_3d::PrimitiveMapT
                                 , rapter::PointContainerT
                                 , rapter::PclCloudT>( int argc, char** argv );
} //...ns rapter
