#ifndef TYPESGL_H
#define TYPESGL_H

#include "types.h"
#include "sampler.h"
#include <QtOpenGL>

namespace InputGen{
namespace Application{
//! Convenience wrapper to call OpenGL command with compile-time vertex definitions
template <typename _Scalar>
struct GLDisplayFunctor{
    static inline void displayVertex(const _Scalar *) {}
};

template <>
inline void
GLDisplayFunctor<double>::displayVertex(const double* data){
    glVertex3dv(data);
}

template <>
inline void
GLDisplayFunctor<float>::displayVertex(const float* data){
    glVertex3fv(data);
}

typedef VisibleSampler<InputGen::Application::Scalar,
                               InputGen::Application::GLDisplayFunctor,
                               InputGen::Application::Primitive>
        Sampler;
}
}


#endif // TYPESGL_H
