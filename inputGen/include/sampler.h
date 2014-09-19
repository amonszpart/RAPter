#ifndef SAMPLER_H
#define SAMPLER_H

#include "primitive.h"

namespace InputGen{


template <typename _Scalar, template <class> class _DisplayFunctor>
class VisibleSampler{
protected:
    typedef _Scalar Scalar;
    typedef _DisplayFunctor<_Scalar> DisplayFunctor;

public:
    virtual void display() const = 0;
    virtual VisibleSampler<_Scalar, _DisplayFunctor>* copy() = 0;

};

//! Basic sampler, sampling the primitive regularly
template <typename _Scalar,
          template <class> class _DisplayFunctor>
struct PrimitiveSampler : public VisibleSampler<
        _Scalar,
        _DisplayFunctor>{
    typedef _Scalar Scalar;
    typedef _DisplayFunctor<_Scalar> DisplayFunctor;

    ///// parameters
    Scalar spacing;


    ///// copy constructor
    /// inline
    PrimitiveSampler(PrimitiveSampler<_Scalar, _DisplayFunctor>* other)
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

    virtual VisibleSampler<_Scalar, _DisplayFunctor>* copy(){
        return new PrimitiveSampler<_Scalar, _DisplayFunctor>(this);
    }
};

////! Basic sampler, sampling the primitive regularly
//template <typename _Scalar,
//          template <class> class _DisplayFunctor>
//struct PonctualSampler : public VisibleSampler<
//        _Scalar,
//        _DisplayFunctor>{
//    typedef _Scalar Scalar;
//    typedef _DisplayFunctor<_Scalar> DisplayFunctor;

//    ///// parameters
//    Scalar angularSpacing;



//    ///// copy constructor
//    /// inline
//    PrimitiveSampler(PrimitiveSampler<_Scalar, _DisplayFunctor>* other)
//        : angularSpacing(angularSpacing->spacing)
//    {}
//    PrimitiveSampler()
//        : angularSpacing(1.)
//    {}

//    ///// processing
//    template <class SampleContainer, class PrimitiveContainer>
//    inline void generateSamples(      SampleContainer&    scontainer,
//                                const PrimitiveContainer& pcontainer);

//    virtual void display() const {
//        _DisplayFunctor dfunctor;
//    }

//    virtual VisibleSampler<_Scalar, _DisplayFunctor>* copy(){
//        return new PonctualSampler<_Scalar, _DisplayFunctor>(this);
//    }
//};


}

#include "impl/sampler.hpp"

#endif // SAMPLER_H
