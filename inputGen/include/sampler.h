#ifndef SAMPLER_H
#define SAMPLER_H

#include "primitive.h"

namespace InputGen{


template <typename _Scalar, template <class> class _DisplayFunctor, class _Primitive>
class VisibleSampler{
protected:
    typedef _Scalar Scalar;
    typedef _Primitive Primitive;
    typedef _DisplayFunctor<_Scalar> DisplayFunctor;

public:
    virtual void display() const = 0;
    virtual VisibleSampler<_Scalar, _DisplayFunctor, _Primitive>* copy() = 0;

protected:

    template <class _SampleContainer, class _Sample>
    inline static
    bool
    addSample(_SampleContainer& scontainer,
              const _Sample& sample,
              const Primitive& prim) {
        if(prim.isInside( sample ) ) // implicit cast should work here
            scontainer.push_back(sample);
            return true;
    }

};

//! Basic sampler, sampling the primitive regularly
template <typename _Scalar,
          template <class> class _DisplayFunctor,
          class _Primitive>
struct PrimitiveSampler : public VisibleSampler<
        _Scalar,
        _DisplayFunctor,
        _Primitive>{
    typedef _Scalar Scalar;
    typedef _Primitive Primitive;
    typedef _DisplayFunctor<_Scalar> DisplayFunctor;
    typedef VisibleSampler< _Scalar, _DisplayFunctor, _Primitive> Base;

    ///// parameters
    Scalar spacing;


    ///// copy constructor
    /// inline
    PrimitiveSampler(PrimitiveSampler<_Scalar, _DisplayFunctor, Primitive>* other)
        : spacing(other->spacing)
    {}
    PrimitiveSampler()
        : spacing(1.)
    {}

    ///// processing
    template <class SampleContainer, class PrimitiveContainer>
    inline void generateSamples(      SampleContainer&    scontainer,
                                const PrimitiveContainer& pcontainer);

    virtual void display() const { }

    virtual VisibleSampler<_Scalar, _DisplayFunctor, _Primitive>* copy(){
        return new PrimitiveSampler<_Scalar, _DisplayFunctor, _Primitive>(this);
    }
};

//! Basic sampler, sampling the primitive regularly
template <typename _Scalar,
          template <class> class _DisplayFunctor,
          class _Primitive>
struct PunctualSampler : public VisibleSampler<
        _Scalar,
        _DisplayFunctor,
        _Primitive>{
    typedef _Scalar Scalar;
    typedef _Primitive Primitive;
    typedef typename _Primitive::vec vec;
    typedef _DisplayFunctor<_Scalar> DisplayFunctor;
    typedef VisibleSampler< _Scalar, _DisplayFunctor, _Primitive> Base;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    ///// parameters
    int nbSamples;  //
    bool occlusion; //! <\brief enable occlusions
    vec pos;



    ///// copy constructor
    /// inline
    PunctualSampler(PunctualSampler<_Scalar, _DisplayFunctor, _Primitive>* other)
        : nbSamples(other->nbSamples),
          occlusion(other->occlusion),
          pos(other->pos)
    {}
    PunctualSampler()
        : nbSamples(1), occlusion(true), pos(vec::Zero())
    {}

    ///// processing
    template <class SampleContainer, class PrimitiveContainer>
    inline void generateSamples(      SampleContainer&    scontainer,
                                const PrimitiveContainer& pcontainer);

    virtual void display() const {
        DisplayFunctor::displayVertex(pos.data());
    }

    virtual VisibleSampler<_Scalar, _DisplayFunctor, _Primitive>* copy(){
        return new PunctualSampler<_Scalar, _DisplayFunctor, _Primitive>(this);
    }
};


}

#include "impl/sampler.hpp"

#endif // SAMPLER_H
