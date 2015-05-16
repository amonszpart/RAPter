#ifndef RAPTER_MERGING_H
#define RAPTER_MERGING_H

#include <iostream>

namespace rapter {

class Merging
{
    public:

        template < class    _PrimitiveContainerT
                 , class    _PointContainerT
                 , typename _Scalar          = float
                 , class    _PointPrimitiveT = typename _PointContainerT::value_type
                 , class    _PrimitiveT      = typename _PrimitiveContainerT::value_type::value_type
                 >
        static inline int mergeCli( int argc, char** argv );

        /*! \brief Greedily assigns points with GID-s that are not in prims to prims that explain them.
        *        Unambiguous assignments go through first, than based on proximity, capped by scale.
        *
        * \tparam _PointPrimitiveDistanceFunctor Has an eval function for a point and all primitives, to calculate the distance from point to primitive. Concept: \ref MyPointPrimitiveDistanceFunctor.
        * \tparam _PointPrimitiveT     Wraps a point, exposing pos() and dir() functions. Concept: \ref rapter::PointPrimitive.
        * \tparam _PrimitiveT          Wraps a primitive, exposing pos() and dir() functions. Concept: \ref rapter::LinePrimitive2.
        * \tparam _PointContainerT     Holds the points. Concept: std::vector< \ref rapter::PointPrimitive >.
        * \tparam _PrimitiveContainerT Holds the primitives grouped by GID in a two level structure. Concept: std::map< int, std::vector< \ref rapter::LinePrimitive2 > >.
        * \tparam _Scalar              Floating point precision of primitives, points, etc. Concept: \ref rapter::PointPrimitive::Scalar.
        * \param[in,out] points        Contains the points, some assigned, some to be assigned to the primitives in prims.
        * \param[in] prims             Contains some primitives tagged with GID and DIR_GID. GID defines the assignment between points and primitives.
        * \param[in] scale             Distance threshold parameter.
        * \param[in] mode              1: re-assign un-ambiguous points (1 adopter); 2: first re-assign unambiguous, then closest, if inside explaining primitive's scale.
        */
        template < class _PointPrimitiveDistanceFunctor
                 , class _PointPrimitiveT
                 , class _PrimitiveT
                 , class _inner_const_iterator
                 , class _PointContainerT
                 , class _PrimitiveContainerT
                 , typename _Scalar >
        static inline int adoptPoints( _PointContainerT &points, _PrimitiveContainerT const& prims, _Scalar const scale, char const mode, int pop_limit );

        /*! \brief Merges adjacent patches that have the same direction ID or are almost parallel.
         *  \tparam _PatchPatchDistanceFunctorT  Concept: \ref rapter::RepresentativeSqrPatchPatchDistanceFunctorT.
         *  \tparam _InnerPrimitiveContainerT    Concept: std::vector<_PrimitiveT>
         *  \param[in] patchPatchDistFunct       Distance functor between two patches, to define adjacency.
         *  \param[in] spatial_threshold         Two extrema should be at least this close to be merged. Concept: \ref MergeParams::spatial_threshold_mult == 3 * scale.
         */
        template < class    _PrimitiveT
                 , class    _PointPrimitiveT
                 //, class    _inner_const_iterator
                 , class    _InnerPrimitiveContainerT
                 , class    _PrimitiveContainerT
                 , class    _PointContainerT
                 , typename _Scalar
                 , class    _PatchPatchDistanceFunctorT
                 , class    _PrimitiveDecideMergeFunctorT
                 , class    _ComparedUidsT
                 >
        static inline int mergeSameDirGids( _PrimitiveContainerT      & out_primitives
                                          , _PointContainerT          & points
                                          , _PrimitiveContainerT /*const&*/ primitives
                                          , _Scalar                       const  scale
                                          , _Scalar                       const  spatial_threshold
                                          , _Scalar                       const  parallel_limit
                                          , _PatchPatchDistanceFunctorT   const& patchPatchDistFunct
                                          , _PrimitiveDecideMergeFunctorT const& primitiveDecideMergeFunct
                                          , bool                                 preserveSmallPatches
                                          , _ComparedUidsT                     & comparedUids
                                          );


}; //...class Merging

} //...namespace RAPTER

#include "rapter/optimization/impl/merging.hpp"

#endif // RAPTER_MERGING_H
