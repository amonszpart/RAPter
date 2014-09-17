#ifndef GF2_PROBLEMSETUP_H
#define GF2_PROBLEMSETUP_H

#include <vector>
#include <map>
#include <set>

#include "qcqpcpp/optProblem.h"  // OptProblem
#include "globfit2/parameters.h" // ProblemSetupParams

namespace GF2
{
    // predecl
    template <typename _Scalar, class _PrimitiveT> struct AbstractPrimitivePrimitiveEnergyFunctor;

    namespace problemSetup
    {
        //! \brief General problem type, most implementations require double, so it is fixed to double.
        typedef qcqpcpp::OptProblem<double> OptProblemT;

        //! \brief              Adds constraints to \p problem so, that each point has at least one line in \p scale radius that is selected.
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

        //! \brief              Adds constraints to \p problem so, that each patch (prims[i] that have the same _PrimitiveT::GID) has at least one member j (prims[i][j]) selected.
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
        //! \param[in,out] problem      The optimization problem to add to.
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
                                           , int                  const  pop_limit
                                           , bool                 const  verbose = false );

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

        //! \brief              Adds unary costs to problem based on point to primitive associations.
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
                                , _Scalar              const  /*scale*/
                                , _Scalar              const  freq_weight );

        //! \brief Nic's version, unfinished!.
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        //! \warning Unfinished
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
    } //... namespace problemSetup

    //! \brief Class to formulate problem into an quadratic optimization problem.
    class ProblemSetup
    {
        public:
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

            /*! \brief                          Step 2. Reads the output from generate and sets up the optimization problem in form of sparse matrices.
             *  \tparam _PrimitiveContainerT    Concept: vector< vector< \ref GF2::LinePrimitive2 > >.
             *  \tparam _PointContainerT        Concept: vector< \ref GF2::PointPrimitive >.
             *  \param[in,out] problem          QCQPcpp::OptProblem (\ref problemSetup::OptProblemT) typed problem to append to.
             *  \param[in] prims                Holds the input primitives to setup the optimization with.
             *  \param[in] points               Holds the input points.
             *  \param[in] constr_mode          See \ref ProblemSetupParams::CONSTR_MODE. \copybrief GF2::ProblemSetupParams::CONSTR_MODE.
             *  \param[in] data_cost_mode       See \ref ProblemSetupParams::DATA_COST_MODE.
             *  \param[in] scale                Scale of input, can be used for data cost.
             *  \param[in] weights              Weights for problemsetup, contains data, pairwise and complexity weight in this order.
             *  \param[in] primPrimDistFunctor  Functor, that calculates the distance between two primitives. Concept: \ref SqrtPrimitivePrimitiveEnergyFunctor.
             *  \param[in] patch_pop_limit      Decides, whether a patch is large or small. Used in \ref problemSetup::largePatchesNeedDirectionConstraint.
             *  \param[in] dir_id_bias          \copydoc ProblemSetupParams::dir_id_bias.
             *  \param[in] verbose              Debug messages display.
             *  \param[in] freq_weight          Multiplies the data cost by freq_weight / DIR_COUNT.
             *  \return                         Outputs EXIT_SUCCESS or the error the OptProblem implementation returns.
             *  \note                           \p points are assumed to be tagged at _PointPrimitiveT::GID with the _PrimitiveT::GID of the \p prims.
             *  \sa \ref problemSetup::largePatchesNeedDirectionConstraint
             */
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
                     , typename ProblemSetupParams<_Scalar>::CONSTR_MODE                  const  constr_mode
                     , typename ProblemSetupParams<_Scalar>::DATA_COST_MODE               const  data_cost_mode
                     , _Scalar                                                            const  scale
                     , Eigen::Matrix<_Scalar,-1,1>                                        const& weights
                     , GF2::AbstractPrimitivePrimitiveEnergyFunctor<_Scalar,_PrimitiveT>* const& primPrimDistFunctor
                     , int                                                                const  patch_pop_limit
                     , _Scalar                                                            const  dir_id_bias
                     , int                                                                const  verbose                = 0
                     , _Scalar                                                            const  freq_weight            = 0.
                     );

    }; //...class ProblemSetup
} //...namespace GF2

#ifndef GF2_INC_PROBLEMSETUP_HPP
#   define GF2_INC_PROBLEMSETUP_HPP
#   include "globfit2/optimization/impl/problemSetup.hpp"
#endif //GF2_INC_PROBLEMSETUP_HPP

#endif // GF2_PROBLEMSETUP_H
