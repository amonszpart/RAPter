#ifndef GF2_SIMPLE_TYPES_H
#define GF2_SIMPLE_TYPES_H

#include <map> // pair
#include <vector>

namespace GF2
{

typedef unsigned int UidT; // UniqueId in containers.hpp
typedef int GidT; // GroupId type
typedef int DidT; // DirectionId type
typedef int LidT; // LinearId type (vector index)
typedef std::pair<GidT,int> GidLid; // uniquely identifies a primitive by first: gid, second: linear id in innerContainer (vector).
typedef std::pair<DidT,int> DidAid; // <DirectionId, AngleId>

typedef float __Scalar; // watch out ,this is usually read from PointPrimitiveT::Scalar elsewhere

struct AnglesT : public std::vector<__Scalar>
{
    typedef __Scalar Scalar;
    AnglesT() {}
    AnglesT( std::vector<Scalar> const& angles ) : AnglesT( angles ) {}
};

} //...ns GF2

#endif // GF2_SIMPLE_TYPES_H
