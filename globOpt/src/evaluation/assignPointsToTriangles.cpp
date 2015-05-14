#include "globfit2/globOpt_types.h"
//#include "globopt/primitives/impl/planePrimitive.hpp"
#include "globopt/evaluation/impl/assignPointsToTriangles.hpp"

namespace globopt
{
    template
    int assignPointsToTriangles< GF2::_3d::PrimitiveVectorT
                                 , GF2::_3d::PrimitiveMapT
                                 , GF2::PointContainerT
                                 , GF2::PclCloudT>( int argc, char** argv );
} //...ns globopt
