#ifndef SAMPLEGENERATOR_H
#define SAMPLEGENERATOR_H

#include "primitive.h"

namespace InputGen{


template <typename _Scalar, template <class> class _DisplayFunctor>
class VisibleSampleGenerator{
protected:
    typedef _Scalar Scalar;
    typedef _DisplayFunctor<_Scalar> DisplayFunctor;

public:
    virtual void display() const = 0;
    virtual VisibleSampleGenerator<_Scalar, _DisplayFunctor>* copy() = 0;

};

/*
 * Generator concept:
 *
\code*/
/*\encode
 *
 * \warning Generators should not have random behaviors, and copied instances
 *          should always generate samples at the same positions.
*/

//! Basic generator, sampling the primitive regularly
template <typename _Scalar,
          template <class> class _DisplayFunctor>
struct PrimitiveSampleGenerator : public VisibleSampleGenerator<
        _Scalar,
        _DisplayFunctor>{
    typedef _Scalar Scalar;
    typedef _DisplayFunctor<_Scalar> DisplayFunctor;

    ///// parameters
    Scalar spacing;


    ///// copy constructor
    /// inline
    PrimitiveSampleGenerator(PrimitiveSampleGenerator<_Scalar, _DisplayFunctor>* other)
        : spacing(other->spacing)
    {}
    PrimitiveSampleGenerator()
        : spacing(1.)
    {}

    ///// processing
    template <class SampleContainer, class PrimitiveContainer>
    inline void generateSamples(      SampleContainer&    scontainer,
                                const PrimitiveContainer& pcontainer);

    virtual void display() const {std::cout <<"display Me !! " << std::endl;}

    virtual VisibleSampleGenerator<_Scalar, _DisplayFunctor>* copy(){
        return new PrimitiveSampleGenerator<_Scalar, _DisplayFunctor>(this);
    }

};



template <typename _Scalar,
          template <class> class T>
template <class SampleContainer, class PrimitiveContainer>
void
PrimitiveSampleGenerator<_Scalar, T>::generateSamples(
              SampleContainer&    scontainer,
        const PrimitiveContainer& pcontainer){
    typedef typename PrimitiveContainer::value_type::vec vec;

    typename PrimitiveContainer::const_iterator it;

    for(it = pcontainer.begin(); it != pcontainer.end(); it++){

        unsigned int nbSampleX = int(std::abs((*it).dim()(0)) / spacing);
        unsigned int nbSampleY = int(std::abs((*it).dim()(1)) / spacing);

        vec midPoint   = (*it).getMidPoint();
        vec tangentVec = (*it).getTangentVector();

        if (nbSampleX == 0 && nbSampleY == 0) // generate a single sample at the primitive midpoint
            scontainer.push_back(midPoint);
        else{

            if (nbSampleX != 0 && nbSampleY == 0){
                scontainer.push_back(midPoint);
                for (unsigned int i = 1; i< nbSampleX/2+1; i++){
                    vec offset = Scalar(i)*spacing*tangentVec;
                    scontainer.push_back(midPoint + offset);
                    scontainer.push_back(midPoint - offset);
                }
                //if ((*it).dim())
            }else if (nbSampleX == 0 && nbSampleY != 0){
                std::cerr << "[Sampler] Unsupported configuration" << std::endl;
            }else{
                std::cerr << "[Sampler] Unsupported configuration" << std::endl;
            }
        }
    }

}


}

#endif // SAMPLEGENERATOR_H
