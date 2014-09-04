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

// direction of the bias. This is updated by sample according to its assignment
vec dir (0., 0., 0.);


switch(biasDirection){
case BIAS_DIRECTION::PRIMITIVE_NORMAL_VECTOR:
{
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
        *darray = dir*bias;
    }

    break;
}
default:
    ;
}

}

#endif // BIASDISPLACEMENT_HPP
