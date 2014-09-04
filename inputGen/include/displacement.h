#ifndef DISPLACEMENT_H
#define DISPLACEMENT_H

#include "primitive.h"
#include <random>
#include <chrono>


namespace InputGen{
    //! Values of this enum are consistent and can be used to set UI
    enum DISPLACEMENT_KERNEL_TYPE{
        DISPLACEMENT_RANDOM_UNIFORM = 0,
        DISPLACEMENT_RANDOM_NORMAL  = 1,
        DISPLACEMENT_BIAS           = 2,
        INVALID_KERNEL              = 3
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

//    template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer>
//    struct RandomDisplacementKernel:
//            public AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>{
//        typedef _Scalar Scalar;
//        virtual void setDistributionRange(Scalar min, Scalar max) = 0;

//        virtual Scalar distributionMin() const = 0;
//        virtual Scalar distributionMax() const = 0;
//    protected:
//        inline
//        RandomDisplacementKernel(const std::string& n,
//                                 DISPLACEMENT_KERNEL_TYPE t)
//            :AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>(n, t) {}

//    };

    template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer, template<typename> class _NumberDistribution>
    struct RandomDisplacementKernel:
            //public RandomDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>{
            public AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>{

    private:
        //! Seed used to generate the samples
        unsigned int _seed;
        std::mt19937 _generator;

    protected:
        typedef _NumberDistribution<_Scalar> NumberDistribution;
        NumberDistribution _distribution;

    public:
        typedef _Scalar Scalar;
        typedef _SampleContainer    SampleContainer;
        typedef _PrimitiveContainer PrimitiveContainer;

        inline RandomDisplacementKernel(const std::string& n,
                                        DISPLACEMENT_KERNEL_TYPE t) :
            AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>(n,t),
            _seed( std::chrono::system_clock::now().time_since_epoch().count()),
            _generator( _seed ),
            _distribution(Scalar(0), Scalar(1))
        {std::cout << "New random kernel " << _seed << std::endl;}

        inline RandomDisplacementKernel(const RandomDisplacementKernel &other) :
            AbstractDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer>(
                other.name,
                other.type),
            _seed( other._seed ),
            _generator( _seed ),
            _distribution( other._distribution )
        {std::cout << "Duplicate random kernel" << _seed << std::endl;}

        virtual void generateDisplacement(
                typename PrimitiveContainer::value_type::vec* darray,
                const SampleContainer& scontainer,
                const PrimitiveContainer& pcontainer);
    };

    /*!
     * Random numbers generated using std::uniform_real_distribution
     */
    template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer>
    struct UniformRandomDisplacementKernel:
            public RandomDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer, std::uniform_real_distribution > {
    protected:
        typedef RandomDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer, std::uniform_real_distribution > Base;
    public:
        inline
        void setDistributionRange(_Scalar min, _Scalar max) {
            Base::_distribution = typename Base::NumberDistribution(min, max);
        }

        inline
        UniformRandomDisplacementKernel()
            : Base( "Random (Uniform)", DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_UNIFORM) {}

        inline _Scalar distributionMin() const { return Base::_distribution.a(); }
        inline _Scalar distributionMax() const { return Base::_distribution.b(); }
    };

    template <typename _Scalar, class _SampleContainer, class _PrimitiveContainer>
    struct NormalRandomDisplacementKernel:
            public RandomDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer, std::normal_distribution > {
    protected:
        typedef RandomDisplacementKernel<_Scalar, _SampleContainer, _PrimitiveContainer, std::normal_distribution > Base;
    public:
        inline
        void setDistributionProperties(_Scalar mean, _Scalar stddev) {
            Base::_distribution = typename Base::NumberDistribution(mean, stddev);
        }

        inline
        NormalRandomDisplacementKernel()
            : Base( "Random (Normal)", DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_NORMAL) {}

        inline _Scalar distributionMean() const { return Base::_distribution.mean(); }
        inline _Scalar distributionStdDev() const { return Base::_distribution.stddev(); }

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
