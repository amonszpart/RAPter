#ifndef SAMPLEGENERATOR_H
#define SAMPLEGENERATOR_H

#include "primitive.h"

namespace InputGen{

/*
template <typename _Scalar>
class SampleGeneratorConcept{

};
*/

template <typename _Scalar>
struct PrimitiveSampleGenerator{
    typedef _Scalar Scalar;

    Scalar spacing;

    template <class SampleContainer, class PrimitiveContainer>
    inline void generateSamples(      SampleContainer&    scontainer,
                                const PrimitiveContainer& pcontainer);

};



template <typename _Scalar>
template <class SampleContainer, class PrimitiveContainer>
void
PrimitiveSampleGenerator<_Scalar>::generateSamples(
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
