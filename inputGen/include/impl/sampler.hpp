#ifndef SAMPLER_HPP
#define SAMPLER_HPP


namespace InputGen{

template <typename _Scalar,
          template <class> class T,
          class _Primitive>
template <class SampleContainer, class PrimitiveContainer>
void
PrimitiveSampler<_Scalar, T, _Primitive>::generateSamples(
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


template<class Derived1, class Derived2>
inline
typename Derived1::Scalar cross2D(const Derived1 & v1, const Derived2 & v2)
{
    return (v1(0)*v2(1)) - (v1(1)*v2(0));
}

template <typename _Scalar,
          template <class> class T,
          class _Primitive>
template <class SampleContainer, class PrimitiveContainer>
void
PunctualSampler<_Scalar, T, _Primitive>::generateSamples(
              SampleContainer&    scontainer,
        const PrimitiveContainer& pcontainer){
    typedef typename PrimitiveContainer::value_type::vec vec;
    typedef typename SampleContainer::value_type Sample;

    typename PrimitiveContainer::const_iterator it;

    _Scalar angleStep = _Scalar(2.) * _Scalar(M_PI) / _Scalar(nbSamples);

    for (unsigned int i = 0; i!= this->nbSamples; i++){
        _Scalar angle = _Scalar(i) * angleStep;
        //direction of the ray from the sample source (this->pos)
        vec d (std::cos(angle), std::sin(angle), _Scalar(0.));

        _Scalar alphaRef = std::numeric_limits<Scalar>::max();
        uint primitiveUID = 0;
        bool found = false;

        // Iterate over all the primitives, compute the intersection with the ray, create a sample on the closest primitive
        for(it = pcontainer.begin(); it != pcontainer.end(); it++){
            // a and b are the end points of the current segment
            const vec a = (*it).coord();
            const vec b = (*it).getEndPoint();

            // we use this technique to compute the intersection
            // http://mathforum.org/library/drmath/view/62814.html

            const _Scalar left  = cross2D(d, a-b);

            if(left != _Scalar(0.)){
                const _Scalar right = cross2D((b-this->pos), a-b);

                // alpha is the distance between the source and the intersection point with the ray
                _Scalar alpha = right / left;

                if (alpha > _Scalar(0.f)){

                    const vec inter = (alpha*d) + this->pos;

                    // now we need to check if we are inside the segment
                    // inter is the position of the intersection, so must be in [0:1] in parametric domain
                    const _Scalar beta  =  (inter - b).norm() / (a-b).norm();
                    if ((inter - b).dot(a-b) >=_Scalar(0.) && beta <= _Scalar(1.)){

                        // now we check for intersection or simply add the segment
                        if (this->occlusion){
                            if (alpha < alphaRef){
                                alphaRef     = alpha;
                                primitiveUID = (*it).uid();
                                found        = true;
                            }
                        }else
                            scontainer.push_back(Sample(inter, alpha*alpha* -d, primitiveUID));
                    }
                }
            }
        }

        if (this->occlusion && found){
            scontainer.push_back(Sample((alphaRef*d) + this->pos, alphaRef*alphaRef*-d, primitiveUID));
        }


            // alpha is the position of the intersection in the ray parametric space, so the distance to the source

            //        unsigned int nbSampleX = int(std::abs((*it).dim()(0)) / spacing);
            //        unsigned int nbSampleY = int(std::abs((*it).dim()(1)) / spacing);

            //        const uint& primitiveUID = (*it).uid();

            //        vec midPoint   = (*it).getMidPoint();
            //        vec tangentVec = (*it).getTangentVector();

            //        if (nbSampleX == 0 && nbSampleY == 0) // generate a single sample at the primitive midpoint
            //            scontainer.push_back(Sample(midPoint, (*it).normal(), primitiveUID));
            //        else{

            //            if (nbSampleX != 0 && nbSampleY == 0){
            //                scontainer.push_back(Sample(midPoint, (*it).normal(), primitiveUID));
            //                for (unsigned int i = 1; i< nbSampleX/2+1; i++){
            //                    vec offset = Scalar(i)*spacing*tangentVec;
            //                    scontainer.push_back(Sample(midPoint + offset, (*it).normal(), primitiveUID));
            //                    scontainer.push_back(Sample(midPoint - offset, (*it).normal(), primitiveUID));
            //                }
            //                //if ((*it).dim())
            //            }else if (nbSampleX == 0 && nbSampleY != 0){
            //                std::cerr << "[Sampler] Unsupported configuration" << std::endl;
            //            }else{
            //                std::cerr << "[Sampler] Unsupported configuration" << std::endl;
            //            }
            //        }

    }

}

}

#endif // SAMPLER_HPP
