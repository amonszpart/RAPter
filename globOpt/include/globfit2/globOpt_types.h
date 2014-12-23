#ifndef GF2_GLOBOPT_TYPES_H
#define GF2_GLOBOPT_TYPES_H

#include <vector>

#include "globfit2/primitives/linePrimitive.h"
#include "globfit2/primitives/planePrimitive.h"
#include "globfit2/primitives/pointPrimitive.h"
#include "globfit2/optimization/energyFunctors.h"

namespace GF2
{
    typedef GF2::PointPrimitive PointPrimitiveT;
    typedef typename PointPrimitiveT::Scalar Scalar;
    typedef std::vector< PointPrimitiveT > PointContainerT;

    namespace _2d
    {
        typedef GF2::LinePrimitive                    PrimitiveT;
        typedef std::vector<PrimitiveT>               InnerPrimitiveContainerT;
        typedef std::vector<InnerPrimitiveContainerT> PrimitiveContainerT;
        typedef MyFinitePrimitiveToFinitePrimitiveCompatFunctor<PrimitiveT/*, MyPointFiniteLineDistanceFunctor*/> MyFiniteLineToFiniteLineCompatFunctor;
    }

    namespace _3d
    {
        typedef GF2::PlanePrimitive                   PrimitiveT;
        typedef std::vector<PrimitiveT>               InnerPrimitiveContainerT;
        typedef std::vector<InnerPrimitiveContainerT> PrimitiveContainerT;
        typedef MyFinitePrimitiveToFinitePrimitiveCompatFunctor<PrimitiveT/*, MyPointFinitePlaneDistanceFunctor*/> MyFinitePlaneToFinitePlaneCompatFunctor;
    }
}

#endif // GF2_GLOBOPT_TYPES_H
