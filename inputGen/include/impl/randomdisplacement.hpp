#ifndef RANDOMDISPLACEMENT_HPP
#define RANDOMDISPLACEMENT_HPP

template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer>
void RandomDisplacementKernel<_Scalar,_SampleContainer,_PrimitiveContainer>::generateDisplacement(
        typename _PrimitiveContainer::value_type::vec* darray,
        const _SampleContainer& scontainer,
        const _PrimitiveContainer& pcontainer){
    typedef typename _PrimitiveContainer::value_type::vec vec;
}

#endif // RANDOMDISPLACEMENT_HPP
