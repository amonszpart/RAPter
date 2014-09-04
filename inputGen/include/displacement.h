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

    template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer>
    struct AbstractDisplacementKernel{
        typedef _SampleContainer    SampleContainer;
        typedef _PrimitiveContainer PrimitiveContainer;

        const std::string name;
        const DISPLACEMENT_KERNEL_TYPE type;

        virtual void generateDisplacement(
                typename PrimitiveContainer::value_type::vec* darray,
                const SampleContainer& scontainer,
                const PrimitiveContainer& pcontainer) = 0;

        protected:
        inline
        AbstractDisplacementKernel(const std::string& n,
                           DISPLACEMENT_KERNEL_TYPE t)
            :name(n), type(t) {}
    };


    template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer>
    struct RandomDisplacementKernel:
            public AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>{
        typedef _Scalar Scalar;
        typedef _SampleContainer    SampleContainer;
        typedef _PrimitiveContainer PrimitiveContainer;

        inline RandomDisplacementKernel() :
            AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>(
                "Random",
                DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM)
        {}

        //! Right now, do nothing
        virtual void generateDisplacement(
                typename PrimitiveContainer::value_type::vec* darray,
                const SampleContainer& scontainer,
                const PrimitiveContainer& pcontainer);
    };


    template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer>
    struct BiasDisplacementKernel:
            public AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>{
        typedef _Scalar Scalar;
        typedef _SampleContainer    SampleContainer;
        typedef _PrimitiveContainer PrimitiveContainer;

        enum BIAS_DIRECTION {
            PRIMITIVE_NORMAL_VECTOR = 0,
            INVALID = 1
        };

        _Scalar bias;
        BIAS_DIRECTION biasDirection;

        inline BiasDisplacementKernel() :
            AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>(
                "Bias",
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
        virtual void generateDisplacement(
                typename PrimitiveContainer::value_type::vec* darray,
                const SampleContainer& scontainer,
                const PrimitiveContainer& pcontainer);
    };


#include "impl/biasdisplacement.hpp"
#include "impl/randomdisplacement.hpp"


}

#endif // DISPLACEMENT_H
