#include "globfit2/optimization/segmentation.h"

#include "globfit2/optimization/segmentation.h" // segmentCli, orientPoints, patchify
#include "globfit2/primitives/pointPrimitive.h"
#include "globfit2/primitives/planePrimitive.h" // segment3D
#include "globfit2/primitives/linePrimitive2.h" // segment2D

int segment3D( int argc, char** argv )
{
    typedef GF2::PlanePrimitive              PrimitiveT;
    typedef std::vector< PrimitiveT >        InnerContainerT;
    typedef std::vector< InnerContainerT >   PrimitiveContainerT; // TODO: map
    typedef GF2::PointPrimitive              PointPrimitiveT;
    typedef std::vector< PointPrimitiveT >   PointContainerT;
    typedef PrimitiveT::Scalar               Scalar;

    return GF2::Segmentation::segmentCli<
                PrimitiveT
                , PrimitiveContainerT
                , PointPrimitiveT
                , PointContainerT
                , Scalar>
            ( argc, argv );
}

int segment2D( int argc, char** argv )
{
    typedef GF2::LinePrimitive2              PrimitiveT;
    typedef std::vector< PrimitiveT >        InnerContainerT;
    typedef std::vector< InnerContainerT > PrimitiveContainerT; // TODO: map
    typedef GF2::PointPrimitive              PointPrimitiveT;
    typedef std::vector< PointPrimitiveT >   PointContainerT;
    typedef PrimitiveT::Scalar               Scalar;

    return GF2::Segmentation::segmentCli<
                PrimitiveT
                , PrimitiveContainerT
                , PointPrimitiveT
                , PointContainerT
                , Scalar>
            ( argc, argv );
}
