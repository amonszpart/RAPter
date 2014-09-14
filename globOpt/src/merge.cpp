#include "globfit2/optimization/merging.h" // mergeCli

#include "optimization/qp/solver.h"

#include "globfit2/primitives/pointPrimitive.h"
#include "globfit2/primitives/planePrimitive.h" // segment3D
#include "globfit2/primitives/linePrimitive2.h" // segment2D


int merge( int argc, char** argv )
{
    return GF2::Merging::mergeCli<GF2::Solver::PrimitiveContainerT,
                                  GF2::Solver::PointContainerT,
                                  GF2::Solver::Scalar>( argc, argv );
}

