#ifndef RAPTER_SIMPLE_TYPES_H
#define RAPTER_SIMPLE_TYPES_H

#include <utility> // pair
#include <vector>

#define RAPTER_MAX_OMP_THREADS (1)

namespace rapter
{
    typedef long                    GidT;   // GroupId type
    typedef unsigned long           UGidT;  // GroupId type
    typedef unsigned long           ULidT;
    typedef GidT                    UidT;   // UniqueId in containers.hpp
    typedef GidT                    DidT;   // DirectionId type // don't ever set to unsigned!!!!
    typedef GidT                    PidT;
    typedef unsigned long           UPidT;
    typedef GidT                    LidT;   // LinearId type (vector index)
    typedef std::pair<GidT,LidT>    GidLid; // uniquely identifies a primitive by first: gid, second: linear id in innerContainer (vector).
    typedef std::pair<DidT,int>     DidAid; // <DirectionId, AngleId>

    typedef float __Scalar; // watch out, this is usually read from PointPrimitiveT::Scalar elsewhere

    const double halfDeg = 0.00872664626; // half degree in radians, used to check generator equality

} //...ns rapter

#endif // RAPTER_SIMPLE_TYPES_H
