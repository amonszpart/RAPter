#ifndef RANDOMDISPLACEMENT_HPP
#define RANDOMDISPLACEMENT_HPP


template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer, template<typename> class NumberDistribution>
void RandomDisplacementKernel<_Scalar,_SampleContainer,_PrimitiveContainer, NumberDistribution>::generateDisplacement(
        typename _PrimitiveContainer::value_type::vec* darray,
        const _SampleContainer& scontainer,
        const _PrimitiveContainer& pcontainer){
    typedef typename _PrimitiveContainer::value_type::vec vec;

    vec dir (0., 0., 0.);

    for (typename SampleContainer::const_iterator it = scontainer.cbegin();
         it != scontainer.cend(); it++, darray++){
        // update direction
        // here we need to iterate over all the primitives to find the one assigned to the point...
        // it is a very slow and unefficient implementation, feel free to change it if you have nothing
        // better to do !
        for (typename PrimitiveContainer::const_iterator itp = pcontainer.cbegin();
             itp != pcontainer.cend(); itp++){
            if ((*itp).uid() == (*it).primitiveId){
                dir = (*itp).normal();
                break;
            }
        }
        // add bias
        *darray = dir*_distribution(_generator);
    }
}

#endif // RANDOMDISPLACEMENT_HPP
