#ifndef TYPES_H
#define TYPES_H

#include "primitive.h"

namespace InputGen{
namespace Application{
typedef double Scalar;
typedef typename InputGen::LinearPrimitive<InputGen::Application::Scalar> Primitive;
}
}

#endif // TYPES_H
