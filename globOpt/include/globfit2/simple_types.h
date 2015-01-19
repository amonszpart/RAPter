#ifndef GF2_SIMPLE_TYPES_H
#define GF2_SIMPLE_TYPES_H

#include <map> // pair
#include <vector>

#define GF2_MAX_OMP_THREADS (2)

namespace GF2
{

typedef long GidT; // GroupId type
typedef unsigned long ULidT;
typedef GidT UidT; // UniqueId in containers.hpp
typedef GidT DidT; // DirectionId type // don't ever set to unsigned!!!!
typedef GidT PidT;
typedef GidT LidT; // LinearId type (vector index)
typedef std::pair<GidT,LidT> GidLid; // uniquely identifies a primitive by first: gid, second: linear id in innerContainer (vector).
typedef std::pair<DidT,int> DidAid; // <DirectionId, AngleId>

typedef float __Scalar; // watch out ,this is usually read from PointPrimitiveT::Scalar elsewhere

struct AnglesT : public std::vector<__Scalar>
{
    typedef __Scalar Scalar;
    typedef std::vector<__Scalar> ParentT;
    AnglesT() {}
    AnglesT( std::vector<Scalar> const& angles ) : ParentT( angles ) {}
};

const double halfDeg = 0.00872664626; // half degree in radians, used to check generator equality

} //...ns GF2

#endif // GF2_SIMPLE_TYPES_H
