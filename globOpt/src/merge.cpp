#include "globfit2/optimization/merging.h"  // mergeCli

#include "globfit2/globOpt_types.h"         // _2d::PrimitiveT, _3d::PrimitiveT
#include "globfit2/util/parse.h"            // find_switch

int merge( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--merge3D") )
    {
        return GF2::Merging::mergeCli< GF2::_3d::PrimitiveContainerT
                                     , GF2::PointContainerT
                                     , GF2::Scalar
                                     , GF2::PointPrimitiveT
                                     , GF2::_3d::PrimitiveT
                                     >( argc, argv );
    }
    else
    {
        return GF2::Merging::mergeCli< GF2::_2d::PrimitiveContainerT
                                     , GF2::PointContainerT
                                     , GF2::Scalar
                                     , GF2::PointPrimitiveT
                                     , GF2::_2d::PrimitiveT
                                     >( argc, argv );
    } //...if find_switch
} //...merge()

