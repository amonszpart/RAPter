#ifndef PROJECT_H
#define PROJECT_H

#include "types.h"
#include "typesGL.h"
#include "displacement.h"

namespace InputGen{
namespace Application{

/*!
 * \brief The Project class
 *
 * This code is a bit messy....
 */

class Project{
private:
    Sampler* _sampler;
public:
    typedef std::vector< InputGen::Application::Primitive > PrimitiveContainer;
    typedef SampleSet SampleContainer;

    PrimitiveContainer primitives;
    SampleContainer    samples;
    typedef std::vector<Primitive::vec,
                        Eigen::aligned_allocator<Primitive::vec> > DisplacementLayer;

    typedef InputGen::AbstractDisplacementKernel<Scalar,
                                 SampleContainer,
                                 PrimitiveContainer> DisplacementKernel;


private:
    //! \brief Structure storing a set of displacement vectors
    struct InternalDisplacementLayer: public DisplacementLayer{

        bool enabled;
        DisplacementKernel * kernel;

        InternalDisplacementLayer( int size,
                                   DisplacementKernel* k,
                                   bool isEnabled = true)
            : enabled(isEnabled), kernel(k)
        { DisplacementLayer::resize(size, Primitive::vec::Zero()); }

        ~InternalDisplacementLayer(){
            // we need smart pointers here, a layer can be copied (and thus destroyed)
            // when a layer is duplicated during layer container resizing
            // \todo Fix that, huge memory leak
            // FIXME Fix that, huge memory leak
            // \warning Need to fix that, huge memory leak
            // delete(kernel);
        }
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
    // Update displacement layer size according to the size of samples
    inline void copySampler(Sampler* tmpSampler)
    {
        delete(_sampler);
        _sampler = tmpSampler->copy();

        for(std::vector<InternalDisplacementLayer>::iterator it = _displ.begin();
            it != _displ.end();
            it++)
            (*it).resize(samples.size());

    }
    inline Sampler* sampler() { return _sampler; }


    //! Noise Layers
    inline void addDisplacementLayer(DisplacementKernel* kernel, bool isEnabled = true)
    { _displ.push_back(InternalDisplacementLayer(samples.size(), kernel, isEnabled)); }

    inline bool isDisplacementLayerEnabled(int layerId) {
        if (layerId >= _displ.size()) return false;
        return _displ[layerId].enabled;
    }

    inline void enableDisplacementLayer(int layerId, bool enabled = true) {
        if (layerId < _displ.size())
            _displ[layerId].enabled = enabled;
    }

    // Ugly, but for now required by DisplacementKernels
    inline Primitive::vec* displacementLayerPtr (int layerId){
        if (layerId >= _displ.size()) return NULL;
        return _displ[layerId].data();
    }

    // Ugly, but for now required by DisplacementKernels
    inline DisplacementKernel* displacementKernel (int layerId){
        if (layerId >= _displ.size()) return NULL;
        return _displ[layerId].kernel;
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
            if ((*it).enabled) displ += (*it)[sampleId];

        return displ;
    }

    inline int nbDisplacementLayers() const { return _displ.size(); }
};
}
}

#endif // PROJECT_H
