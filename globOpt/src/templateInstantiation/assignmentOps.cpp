#include "globfit2/globOpt_types.h"
#include "globopt/processing/impl/assignmentOps.hpp"

namespace globopt
{
    template
    int approxUnassignedWPlanes< GF2::_3d::PrimitiveVectorT
                           , GF2::_3d::PrimitiveMapT
                           , GF2::PointContainerT
                           , GF2::PclCloudT>( int argc, char** argv );
}
