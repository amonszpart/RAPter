#include "globfit2/comparison/impl/globFit.hpp"
//#include "rapter/processing/impl/subsamplePrimitives.hpp"
#include "globfit2/globOpt_types.h"

namespace rapter
{
    template
    int toGlobFit< rapter::PclCloudT
             , rapter::_3d::PrimitiveMapT
             , rapter::_3d::PrimitiveVectorT
             , rapter::PointContainerT>( int argc, char **argv );

    template
    int fromGlobFit< rapter::PclCloudT
                   , rapter::_3d::PrimitiveMapT
                   , rapter::_3d::PrimitiveVectorT
                   , rapter::PointContainerT>( int argc, char **argv );

//    template
//    int subsamplePrimitives< GF2::_3d::PrimitiveVectorT
//                           , GF2::_3d::PrimitiveMapT
//                           , GF2::PointContainerT
//                           , GF2::PclCloudT>( int argc, char** argv );
}
