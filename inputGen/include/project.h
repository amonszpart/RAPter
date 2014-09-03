#ifndef PROJECT_H
#define PROJECT_H

#include "types.h"
#include "typesGL.h"

namespace InputGen{
namespace Application{
class Project{
public:
    std::vector< InputGen::Application::Primitive > primitives;
    InputGen::Application::PointSet samples;

protected:
    Sampler* _sampler;

public:
    inline Project() : _sampler(NULL) {}

    inline void clear () {
        primitives.clear();
        samples.clear();

    }

    //! Store a copy of tmpSampler
    inline void copySampler(Sampler* tmpSampler) { _sampler = tmpSampler->copy(); }
    inline Sampler* sampler() { return _sampler; }
};
}
}

#endif // PROJECT_H
