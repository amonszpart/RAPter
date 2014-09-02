#ifndef TYPES_H
#define TYPES_H

#include "primitive.h"

#include <vector>
#include "eigen3/Eigen/StdVector"

namespace InputGen{
namespace Application{
typedef double Scalar;
typedef typename InputGen::LinearPrimitive<InputGen::Application::Scalar> Primitive;

typedef std::vector<Primitive::vec,
                    Eigen::aligned_allocator<Primitive::vec> > PointSet;
}
}

#endif // TYPES_H
