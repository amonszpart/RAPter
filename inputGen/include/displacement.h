#ifndef DISPLACEMENT_H
#define DISPLACEMENT_H

#include "primitive.h"


namespace InputGen{
    //! Values of this enum are consistent and can be used to set UI
    enum DISPLACEMENT_KERNEL_TYPE{
        DISPLACEMENT_RANDOM = 0,
        DISPLACEMENT_BIAS   = 1,
        INVALID_KERNEL      = 2
    };

    template <typename _Scalar>
    struct DisplacementKernel{
        const std::string name;
        const DISPLACEMENT_KERNEL_TYPE type;

        protected:
        inline
        DisplacementKernel(const std::string& n,
                           DISPLACEMENT_KERNEL_TYPE t)
            :name(n), type(t) {}
    };


    template <typename _Scalar>
    struct RandomDisplacementKernel: public DisplacementKernel<_Scalar>{
        typedef _Scalar Scalar;

        inline RandomDisplacementKernel() :
            DisplacementKernel<_Scalar>("Random",
                                        DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM)
        {}

        //! Right now, do nothing
        template <class SampleContainer, class PrimitiveContainer>
        inline void generateDisplacement( typename PrimitiveContainer::value_type::vec* darray,
                                          const SampleContainer& scontainer,
                                          const PrimitiveContainer& pcontainer) {}
    };


    template <typename _Scalar>
    struct BiasDisplacementKernel: public DisplacementKernel<_Scalar>{
        typedef _Scalar Scalar;

        enum BIAS_DIRECTION {
            PRIMITIVE_NORMAL_VECTOR = 0,
            INVALID = 1
        };

        _Scalar bias;
        BIAS_DIRECTION biasDirection;

        inline BiasDisplacementKernel() :
            DisplacementKernel<_Scalar>("Bias",
                                        DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_BIAS),
            bias(0),
            biasDirection(PRIMITIVE_NORMAL_VECTOR)
        {}


        /*!
         * Todo: would require here a class to access and edit elements of a vector
         * without changing its size
         *
         * Here, darray is ensured to have the same number of elements that scontainer
         */
        template <class SampleContainer, class PrimitiveContainer>
        inline void generateDisplacement( typename PrimitiveContainer::value_type::vec* darray,
                                          const SampleContainer& scontainer,
                                          const PrimitiveContainer& pcontainer);
    };

template <typename _Scalar>
template <class SampleContainer, class PrimitiveContainer>
void BiasDisplacementKernel<_Scalar>::generateDisplacement(
        typename PrimitiveContainer::value_type::vec* darray,
        const SampleContainer& scontainer,
        const PrimitiveContainer& pcontainer){
    typedef typename PrimitiveContainer::value_type::vec vec;

    if (this->biasDirection >= BIAS_DIRECTION::INVALID)
        return;

    // direction of the bias. This is updated by sample according to its assignment
    // This is a very slow and unefficient manner, feel free to change it if you have nothing
    // better to do !
    vec dir (0., 0., 0.);


    for (typename SampleContainer::const_iterator it = scontainer.cbegin();
         it != scontainer.cend(); it++, darray++){
        // update direction
        switch(biasDirection){
        case BIAS_DIRECTION::PRIMITIVE_NORMAL_VECTOR:
        {
            // here we need to iterate over all the primitives to find the right one...
            // unefficient you said ??
            for (typename PrimitiveContainer::const_iterator itp = pcontainer.cbegin();
                 itp != pcontainer.cend(); itp++){
                if ((*itp).uid() == (*it).primitiveId){
                    dir = (*itp).normal();
                    break;
                }
            }
            break;
        }
        }

        // add bias
        *darray = dir*bias;
    }

}

}

#endif // DISPLACEMENT_H
