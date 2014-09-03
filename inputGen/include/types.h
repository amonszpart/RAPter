#ifndef TYPES_H
#define TYPES_H

#include "primitive.h"

#include <vector>
#include "eigen3/Eigen/StdVector"

namespace InputGen{
namespace Application{
typedef double Scalar;
typedef typename InputGen::LinearPrimitive<InputGen::Application::Scalar> Primitive;


//! \brief sample storing a position and its assignment
struct Sample: public Primitive::vec{
    typedef Primitive::vec Base;
    typedef Base::Scalar Scalar;
    int primitiveId;

    inline Sample(): Base(), primitiveId(-1) {}
    template <class Derived>
    inline Sample(const Derived&v): Base(v), primitiveId(-1) {}
    inline Sample(const Scalar&x,
                  const Scalar&y,
                  const Scalar&z): Base(x,y,z), primitiveId(-1) {}
};

typedef std::vector<Sample,
                    Eigen::aligned_allocator<Sample> > PointSet;
}
}

#endif // TYPES_H
