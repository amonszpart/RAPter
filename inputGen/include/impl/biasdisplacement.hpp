#ifndef BIASDISPLACEMENT_HPP
#define BIASDISPLACEMENT_HPP


template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer>
void BiasDisplacementKernel<_Scalar,_SampleContainer,_PrimitiveContainer>::generateDisplacement(
        typename _PrimitiveContainer::value_type::vec* darray,
        const _SampleContainer& scontainer,
        const _PrimitiveContainer& pcontainer){
    typedef typename _PrimitiveContainer::value_type::vec vec;

    if (this->biasDirection >= BIAS_DIRECTION::INVALID)
        return;

    switch(biasDirection){
    case BIAS_DIRECTION::PRIMITIVE_NORMAL_VECTOR:
    {
        for (typename SampleContainer::const_iterator it = scontainer.cbegin();
             it != scontainer.cend(); it++, darray++){
            // add bias
            *darray = (*it).normal*bias;
        }

        break;
    }
    default:
        ;
    }

}

#endif // BIASDISPLACEMENT_HPP
