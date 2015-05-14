#include "globfit2/globOpt_types.h"
#include "globopt/processing/impl/subsamplePrimitives.hpp"

namespace globopt
{
    template
    int subsamplePrimitives< GF2::_3d::PrimitiveVectorT
                           , GF2::_3d::PrimitiveMapT
                           , GF2::PointContainerT
                           , GF2::PclCloudT>
                           ( int argc, char** argv );

    template
    int subsamplePrimitives< GF2::_2d::PrimitiveVectorT
                           , GF2::_2d::PrimitiveMapT
                           , GF2::PointContainerT
                           , GF2::PclCloudT>
                           ( int argc, char** argv );
} //...ns globopt
