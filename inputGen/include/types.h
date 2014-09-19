#ifndef TYPES_H
#define TYPES_H

#include "primitive.h"

#include <vector>
#include "Eigen/StdVector"

namespace InputGen{
namespace Application{
typedef double Scalar;
typedef typename InputGen::LinearPrimitive<InputGen::Application::Scalar> Primitive;


//! \brief sample storing a position and its assignment
struct Sample: public Primitive::vec{
    typedef Primitive::vec Base;
    typedef Base::Scalar Scalar;
    int primitiveId;
    Primitive::vec normal; // normal is mainly used to compute displacement

    inline Sample(int id = -1): Base(), primitiveId(id), normal(Base::Zero()) {}
    template <class Derived1, class Derived2>
    inline Sample(const Derived1&v, const Derived2&n, int id = -1): Base(v), primitiveId(id), normal(n) {}
    inline Sample(const Scalar&x,
                  const Scalar&y,
                  const Scalar&z,
                  int id = -1): Base(x,y,z), primitiveId(id), normal(Base::Zero()) {}
    inline Sample(const Scalar&x,
                  const Scalar&y,
                  const Scalar&z,
                  const Scalar&nx,
                  const Scalar&ny,
                  const Scalar&nz,
                  int id = -1): Base(x,y,z), primitiveId(id), normal(nx, ny, nz) {}
};

typedef std::vector<Sample,
                    Eigen::aligned_allocator<Sample> > SampleSet;
}

template <typename Scalar>
struct MergeParam{
    bool useAngular;
    Scalar angleRef; // stored in radian
    bool usePeriodicAngles;
};
}


#endif // TYPES_H
