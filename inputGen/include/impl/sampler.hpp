#ifndef SAMPLER_HPP
#define SAMPLER_HPP


namespace InputGen{

template <typename _Scalar,
          template <class> class T>
template <class SampleContainer, class PrimitiveContainer>
void
PrimitiveSampler<_Scalar, T>::generateSamples(
              SampleContainer&    scontainer,
        const PrimitiveContainer& pcontainer){
    typedef typename PrimitiveContainer::value_type::vec vec;
    typedef typename SampleContainer::value_type Sample;

    typename PrimitiveContainer::const_iterator it;

    // Iterate over all the primitives and generate samples using uniform spacing
    for(it = pcontainer.begin(); it != pcontainer.end(); it++){

        unsigned int nbSampleX = int(std::abs((*it).dim()(0)) / spacing);
        unsigned int nbSampleY = int(std::abs((*it).dim()(1)) / spacing);

        const uint& primitiveUID = (*it).uid();

        vec midPoint   = (*it).getMidPoint();
        vec tangentVec = (*it).getTangentVector();

        if (nbSampleX == 0 && nbSampleY == 0) // generate a single sample at the primitive midpoint
            scontainer.push_back(Sample(midPoint, (*it).normal(), primitiveUID));
        else{

            if (nbSampleX != 0 && nbSampleY == 0){
                scontainer.push_back(Sample(midPoint, (*it).normal(), primitiveUID));
                for (unsigned int i = 1; i< nbSampleX/2+1; i++){
                    vec offset = Scalar(i)*spacing*tangentVec;
                    scontainer.push_back(Sample(midPoint + offset, (*it).normal(), primitiveUID));
                    scontainer.push_back(Sample(midPoint - offset, (*it).normal(), primitiveUID));
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

#endif // SAMPLER_HPP
