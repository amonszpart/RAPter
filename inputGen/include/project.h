#ifndef PROJECT_H
#define PROJECT_H

#include "types.h"
#include "typesGL.h"

namespace InputGen{
namespace Application{

class Project{
private:
    Sampler* _sampler;
public:
    std::vector< InputGen::Application::Primitive > primitives;
    InputGen::Application::SampleSet samples;
    typedef std::vector<Primitive::vec,
                        Eigen::aligned_allocator<Primitive::vec> > DisplacementLayer;


private:
    //! \brief Structure storing a set of displacement vectors
    struct InternalDisplacementLayer: public DisplacementLayer{
        bool enabled;

        InternalDisplacementLayer( int size,
                            bool isEnabled = true)
            : enabled(isEnabled)
        { DisplacementLayer::resize(size, Primitive::vec::Zero()); }
    };

    const DisplacementLayer::iterator _nullDisplacementIterator;
    std::vector<InternalDisplacementLayer> _displ;

public:
    inline Project() :
        _sampler(NULL),
        _nullDisplacementIterator(DisplacementLayer().begin()) {}

    inline void clear () {
        primitives.clear();
        samples.clear();
        _displ.clear();
        delete(_sampler);
    }

    //! Store a copy of tmpSampler
    inline void copySampler(Sampler* tmpSampler) { _sampler = tmpSampler->copy(); }
    inline Sampler* sampler() { return _sampler; }


    //! Noise Layers
    inline void addDisplacementLayer(bool isEnabled)
    { _displ.push_back(InternalDisplacementLayer(samples.size(), isEnabled)); }

    inline bool isDisplacementLayerEnabled(int layerId) {
        if (layerId >= _displ.size()) return false;
        return _displ[layerId].enabled;
    }

    inline DisplacementLayer::iterator beginDisplacementLayer(int layerId) {
        if (layerId >= _displ.size()) return _nullDisplacementIterator;
        return _displ[layerId].begin();
    }

    inline DisplacementLayer::iterator endDisplacementLayer(int layerId) {
        if (layerId >= _displ.size()) return _nullDisplacementIterator;
        return _displ[layerId].end();
    }

    inline DisplacementLayer::const_iterator cbeginDisplacementLayer(int layerId) const {
        if (layerId >= _displ.size()) return _nullDisplacementIterator;
        return _displ[layerId].cbegin();
    }

    inline DisplacementLayer::const_iterator cendDisplacementLayer(int layerId) const {
        if (layerId >= _displ.size()) return _nullDisplacementIterator;
        return _displ[layerId].cend();
    }

    // we do not use an iterator here to emphasis that this value is computed
    inline Primitive::vec computeTotalDisplacement(int sampleId) const {
        Primitive::vec displ (Primitive::vec::Zero());
        for(std::vector<InternalDisplacementLayer>::const_iterator it = _displ.begin();
            it != _displ.end();
            it++)
            displ += (*it)[sampleId];

        return displ;
    }
};
}
}

#endif // PROJECT_H
