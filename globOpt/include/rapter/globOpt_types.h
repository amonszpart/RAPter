#ifndef RAPTER_GLOBOPT_TYPES_H
#define RAPTER_GLOBOPT_TYPES_H

#include <vector>

#include "rapter/simple_types.h"
#include "rapter/primitives/linePrimitive.h"
#include "rapter/primitives/planePrimitive.h"
#include "rapter/primitives/pointPrimitive.h"
#include "rapter/optimization/energyFunctors.h"
#include "rapter/util/containers.hpp"

namespace rapter
{
    typedef rapter::PointPrimitive PointPrimitiveT;
    //typedef typename PointPrimitiveT::Scalar Scalar;
    typedef __Scalar Scalar;
    //typedef PrimitiveVectorT< PointPrimitiveT > PointContainerT;
    typedef PointPrimitiveVector PointContainerT;

    namespace _2d
    {
        typedef rapter::LinePrimitive                    PrimitiveT;
        typedef std::vector<PrimitiveT>               InnerPrimitiveContainerT;
        //typedef std::vector<InnerPrimitiveContainerT> PrimitiveContainerT;
        typedef PrimitiveVectorT<PrimitiveT>          PrimitiveContainerT;
        typedef PrimitiveVectorT<PrimitiveT>          PrimitiveVectorT;
        typedef containers::PrimitiveContainer<PrimitiveT>          PrimitiveMapT;
        typedef MyFinitePrimitiveToFinitePrimitiveCompatFunctor<PrimitiveT/*, MyPointFiniteLineDistanceFunctor*/> MyFiniteLineToFiniteLineCompatFunctor;
    }

    namespace _3d
    {
        typedef rapter::PlanePrimitive                   PrimitiveT;
        typedef std::vector<PrimitiveT>               InnerPrimitiveContainerT;
        //typedef std::vector<InnerPrimitiveContainerT> PrimitiveContainerT;
        typedef PrimitiveVectorT<PrimitiveT>          PrimitiveContainerT;
        typedef PrimitiveVectorT<PrimitiveT>          PrimitiveVectorT;
        typedef containers::PrimitiveContainer<PrimitiveT>          PrimitiveMapT;
        typedef MyFinitePrimitiveToFinitePrimitiveCompatFunctor<PrimitiveT/*, MyPointFinitePlaneDistanceFunctor*/> MyFinitePlaneToFinitePlaneCompatFunctor;
    }
}

#endif // RAPTER_GLOBOPT_TYPES_H
