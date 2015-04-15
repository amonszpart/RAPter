#include "globfit2/globOpt_types.h"
#include "globopt/processing/impl/assignmentOps.hpp"
#include "globfit2/primitives/planePrimitive.h"

namespace globopt
{
    template
    int approxUnassignedWPlanes< GF2::_3d::PrimitiveVectorT
                           , GF2::_3d::PrimitiveMapT
                           , GF2::PointContainerT
                           , GF2::PclCloudT>( int argc, char** argv );
}
