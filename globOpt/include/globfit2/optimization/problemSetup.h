#ifndef GF2_PROBLEMSETUP_H
#define GF2_PROBLEMSETUP_H

#include <vector>
#include <map>
#include <set>

#include "qcqpcpp/optProblem.h"  // OptProblem
//#include "globfit2/parameters.h" // ProblemSetupParams

namespace GF2
{
    // predecl
    template <typename _Scalar, class _PrimitiveT> struct AbstractPrimitivePrimitiveEnergyFunctor;

    namespace problemSetup
    {
        typedef qcqpcpp::OptProblem<double> OptProblemT;

        //! \brief Adds unary costs to problem based on point to primitive associations.
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        associationBasedDataCost( _OptProblemT              & problem
                                , _PrimitiveContainerT const& prims
                                , _PointContainerT     const& points
                                , _AssocT              const& lids_varids
                                , _WeightsT            const& weights
                                , _Scalar              const  /*scale*/ );

        //! \brief
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        everyPatchNeedsDirectionConstraint( _OptProblemT              & problem
                                          , _PrimitiveContainerT const& prims
                                          , _PointContainerT     const& /*points*/
                                          , _AssocT              const& lids_varids
                                          , _WeightsT            const& /*weights*/
                                          , _Scalar              const  /*scale*/ );

        //! \brief                      Hybrid mode, large patches need one direction, points in small patches need to be assigned to one.
        //! \tparam _PrimitiveT         Concept: _PrimitiveContainerT::value_type::value_type
        //! \tparam _PointPrimitiveT    Concept: _PointContainerT::value_type
        //! \tparam _AssocT             Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        //! \param[in/out] problem      The optimization problem to add to.
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT
                 , class _PointPrimitiveT
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        largePatchesNeedDirectionConstraint( _OptProblemT              & problem
                                           , _PrimitiveContainerT const& prims
                                           , _PointContainerT     const& points
                                           , _AssocT              const& lids_varids
                                           , _WeightsT            const& /*weights*/
                                           , _Scalar              const  scale
                                           , int                  const  pop_limit );


        //! \brief Nic's version
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        instanceBasedDataCost( _OptProblemT              & problem
                                , _PrimitiveContainerT const& prims
                                , _PointContainerT     const& points
                                , _AssocT              const& lids_varids
                                , _WeightsT            const& weights
                                , _Scalar              const  /*scale*/ );


        //! \brief Adds unary costs to problem based on scale wide band assocation.
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        bandBasedDataCost( _OptProblemT              & problem
                         , _PrimitiveContainerT const& prims
                         , _PointContainerT     const& points
                         , _AssocT              const& lids_varids
                         , _WeightsT            const& weights
                         , _Scalar              const  scale );

        //! \brief
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        everyPointNeedsPatchConstraint( _OptProblemT              & problem
                                      , _PrimitiveContainerT const& prims
                                      , _PointContainerT     const& points
                                      , _AssocT              const& lids_varids
                                      , _WeightsT            const& /*weights*/
                                      , _Scalar              const  scale );

        //! \brief doxytest
        static int dummy() {}
    } //... namespace problemSetup

    //! \brief Class to formulate problem into an quadratic optimization problem.
    class ProblemSetup
    {
        public:
            //! \brief Options to calculate data cost of primitive by.
            enum DATA_COST_MODE { ASSOC_BASED    = 0 //!< \brief \sa \ref problemSetup::associationBasedDataCost
                                , BAND_BASED     = 1 //!< \brief \sa \ref problemSetup::bandBasedDataCost
                                , INSTANCE_BASED = 2 //!< \brief \sa \ref problemSetup::instanceBasedDataCost
                                };
            enum CONSTR_MODE    { PATCH_WISE = 0, POINT_WISE = 1, HYBRID = 2 };

            //! \brief                          Step 2. Reads the output from generate and sets up the optimization problem in form of sparse matrices.
            //! \tparam _PrimitiveContainerT    Concept: vector< vector< \ref GF2::LinePrimitive2 > >
            //! \tparam _PointContainerT        Concept: vector< \ref GF2::PointPrimitive >
            //! \param argc                     Number of CLI arguments.
            //! \param argv                     Vector of CLI arguments.
            //! \return                         Outputs EXIT_SUCCESS or the error the OptProblem implementation returns.
            //! \sa \ref problemSetup::largePatchesNeedDirectionConstraint
            template <
                       class _PrimitiveContainerT
                     , class _PointContainerT
                     , class _PrimitiveT          = typename _PrimitiveContainerT::value_type::value_type
                     , class _PointPrimitiveT     = typename _PointContainerT::value_type
                     >
            static inline int formulateCli( int argc, char** argv );

            //! \brief                          Step 2. Reads the output from generate and sets up the optimization problem in form of sparse matrices.
            //! \tparam _PrimitiveContainerT    Concept: vector< vector< \ref GF2::LinePrimitive2 > >.
            //! \tparam _PointContainerT        Concept: vector< \ref GF2::PointPrimitive >.
            //! \param  prims                   Holds the input primitives to setup the optimization with.
            //! \param  points                  Holds the input points.
            //! \param  primPrimDistFunctor     Functor, that calculates the distance between two primitives. Concept: \ref SqrtPrimitivePrimitiveEnergyFunctor.
            //! \return                         Outputs EXIT_SUCCESS or the error the OptProblem implementation returns.
            //! \note                           \p points are assumed to be tagged at _PointPrimitiveT::GID with the _PrimitiveT::GID of the \p prims.
            //! \sa \ref problemSetup::largePatchesNeedDirectionConstraint
            template < class _PointPrimitiveDistanceFunctor
                     , class _PrimitiveContainerT
                     , class _PointContainerT
                     , class _PrimitiveT        = typename _PrimitiveContainerT::value_type::value_type
                     , class _PointPrimitiveT   = typename _PointContainerT::value_type
                     , typename _Scalar         = typename _PointPrimitiveT::Scalar
                     > static inline int
            formulate( problemSetup::OptProblemT                                               & problem
                     , _PrimitiveContainerT                                               const& prims
                     , _PointContainerT                                                   const& points
                     //, std::vector<std::pair<int,int> >                                   const& points_primitives
                     , CONSTR_MODE                                                        const  constr_mode
                     , DATA_COST_MODE                                                     const  data_cost_mode
                     , _Scalar                                                            const  scale
                     , Eigen::Matrix<_Scalar,-1,1>                                        const& weights
                     , GF2::AbstractPrimitivePrimitiveEnergyFunctor<_Scalar,_PrimitiveT>* const& primPrimDistFunctor
                     , int                                                                const  patch_pop_limit
                     , int                                                                const  verbose = 0 );

    }; //...class ProblemSetup
} //...namespace GF2

#ifndef GF2_INC_PROBLEMSETUP_HPP
#   define GF2_INC_PROBLEMSETUP_HPP
#   include "globfit2/optimization/impl/problemSetup.hpp"
#endif //GF2_INC_PROBLEMSETUP_HPP

#endif // GF2_PROBLEMSETUP_H
