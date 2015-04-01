#include "globfit2/comparison/impl/globFit.hpp"
#include "globfit2/globOpt_types.h"

namespace globopt
{
    template
    int toGlobFit< GF2::PclCloudT
             , GF2::_3d::PrimitiveMapT
             , GF2::_3d::PrimitiveVectorT
             , GF2::PointContainerT>( int argc, char **argv );

    template
    int fromGlobFit< GF2::PclCloudT
                     , GF2::_3d::PrimitiveMapT
                     , GF2::_3d::PrimitiveVectorT
                     , GF2::PointContainerT>( int argc, char **argv );
}
