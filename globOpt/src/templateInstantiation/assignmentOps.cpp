#include "rapter/typedefs.h"
#include "rapter/util/pclUtil.h"                    // PclCloudT
#include "rapter/processing/impl/assignmentOps.hpp"
#include "rapter/primitives/planePrimitive.h"
#include "rapter/primitives/impl/planePrimitive.hpp"

namespace rapter
{
    template
    int approxUnassignedWPlanes< rapter::_3d::PrimitiveVectorT
                               , rapter::_3d::PrimitiveMapT
                               , rapter::PointContainerT
                               , rapter::PclCloudT>( int argc, char** argv );
}
