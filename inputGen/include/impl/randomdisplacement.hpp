#ifndef RANDOMDISPLACEMENT_HPP
#define RANDOMDISPLACEMENT_HPP


template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer, template<typename> class NumberDistribution>
void RandomDisplacementKernel<_Scalar,_SampleContainer,_PrimitiveContainer, NumberDistribution>::generateDisplacement(
        typename _PrimitiveContainer::value_type::vec* darray,
        const _SampleContainer& scontainer,
        const _PrimitiveContainer& pcontainer){
    typedef typename _PrimitiveContainer::value_type::vec vec;

    for (typename SampleContainer::const_iterator it = scontainer.cbegin();
         it != scontainer.cend(); it++, darray++){
        *darray = (*it).normal*_distribution(_generator);
    }
}

#endif // RANDOMDISPLACEMENT_HPP
