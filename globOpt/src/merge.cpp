#include "globfit2/optimization/merging.h" // mergeCli

//#include "optimization/qp/solver.h"

//#include "globfit2/primitives/pointPrimitive.h"
//#include "globfit2/primitives/planePrimitive.h" // segment3D
//#include "globfit2/primitives/linePrimitive2.h" // segment2D

#include "globfit2/globOpt_types.h" // _2d::PrimitiveT, _3d::PrimitiveT
#include "globfit2/util/parse.h" // find_switch

int merge( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--merge3D") )
    {
        std::cerr << "[" << __func__ << "]: " << "merge3D does not compile yet..." << std::endl;
        return EXIT_FAILURE;
        // problem at mergingFunctors line 41, when calling PlanePrimitive::normal
//        return GF2::Merging::mergeCli< GF2::_3d::PrimitiveContainerT
//                                     , GF2::PointContainerT
//                                     , GF2::Scalar
//                                     , GF2::PointPrimitiveT
//                                     , GF2::_3d::PrimitiveT
//                                     >( argc, argv );
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

