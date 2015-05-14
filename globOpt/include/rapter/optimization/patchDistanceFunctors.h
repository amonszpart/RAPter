#ifndef RAPTER_PATCHDISTANCEFUNCTORS_H
#define RAPTER_PATCHDISTANCEFUNCTORS_H

#include "rapter/my_types.h" // angleInRad

namespace rapter
{

#if RAPTER_WITH_FULL_LINKAGE
//! \brief Spatial patch-patch distance is the maximum distance between any two members (farthest two points).
template <typename _Scalar>
struct SpatialPatchPatchMaxDistanceFunctorT
{
        // spatial distanceToPatch
        // Concept: PointContainerT = std::vector<PointPrimitive>, PatchT = std::vector<int>, PointT = PointPrimitive
        template <class _PointT /* = typename _PointContainerT::value_type*/, class _PatchAT, class _PatchBT, class _PointContainerT >
        static inline _Scalar eval( _PatchAT const& patch0, _PatchBT const& patch1, _PointContainerT const& points, _Scalar const* cut_off )
        {
            typedef typename _PointT::VectorType VectorType; // concept: Eigen::Vector3f

            _Scalar max_distance = _Scalar(0);
            for ( int pid_id0 = 0; pid_id0 != patch0.size(); ++pid_id0 )
            {
                const int pid0 = patch0[pid_id0].first;
                for ( int pid_id1 = 0; pid_id1 != patch1.size(); ++pid_id1 )
                {
                    const int pid1 = patch1[ pid_id1 ].first;
                    _Scalar dist  = ((VectorType)points.at(pid1) - (VectorType)points.at(pid0)).norm();
                    if ( dist > max_distance )
                        max_distance = dist;

                    // early stop
                    if ( cut_off && (max_distance > *cut_off ) )
                        return max_distance;
                }
            }

            return max_distance;
        }

        static std::string toString() { return "SpatialPatchPatchMaxDistanceFunctorT"; }
}; //...struct SpatialPatchPatchMaxDistanceFunctorT

//! \brief Spatial patch-patch distance is the minimum distance between any two members (closest two points).
template <typename _Scalar>
struct SpatialPatchPatchMinDistanceFunctorT
{
        // Concept: PointContainerT = std::vector<PointPrimitive>, PatchT = std::vector<int>, PointT = PointPrimitive
        template <class _PointT /*= typename _PointContainerT::value_type*/, class _PatchAT, class _PatchBT, class _PointContainerT>
        static inline _Scalar eval( _PatchAT const& patch0, _PatchBT const& patch1, _PointContainerT const& points, _Scalar const* /*cut_off*/ )
        {
            typedef typename _PointT::VectorType VectorType; // concept: Eigen::Vector3f

            _Scalar min_distance = _Scalar( std::numeric_limits<_Scalar>::max() );
            for ( int pid_id0 = 0; pid_id0 != patch0.size(); ++pid_id0 )
            {
                const int pid0 = patch0[pid_id0].first;
                for ( int pid_id1 = 0; pid_id1 != patch1.size(); ++pid_id1 )
                {
                    const int pid1 = patch1[ pid_id1 ].first;
                    _Scalar dist  = ((VectorType)points.at(pid1) - (VectorType)points.at(pid0)).norm();
                    if ( dist < min_distance )
                        min_distance = dist;
                }
            }
            return min_distance;
        } // ... SpatialPointPatchMinDistanceFunctorT::eval()

        static std::string toString() { return "SpatialPatchPatchMinDistanceFunctorT"; }
};  //...struct SpatialPatchPatchMinDistanceFunctorT
#endif

//! \brief Spatial patch-patch distance is the distance of the two patch positions.
template <typename _Scalar>
struct SpatialPatchPatchSingleDistanceFunctorT
{
        // spatial distanceToPatch
        // Concept: PointT = PointPrimitive, PatchAT, PatchBT = std::vector<int>
        template <class _PointT, class _PatchAT, class _PatchBT, class _PointContainerT>
        static inline _Scalar eval( _PatchAT const& patch0, _PatchBT const& patch1, _PointContainerT const& points, _Scalar const* /*cut_off*/  )
        {
           return (patch0.template pos() - patch1.template pos()).norm();
        }

        static std::string toString() { return "SpatialPatchPatchSingleDistanceFunctorT"; }
}; //...struct SpatialPatchPatchMaxDistanceFunctorT

//_________________________________________________________________________________________________________________________________________________________________________________

//! \brief Calculates and merges distance in spatial and angular domain between two patches.
template <typename _Scalar, class _SpatialPatchPatchDistanceFunctorT>
struct AbstractPatchPatchDistanceFunctorT
{
    inline AbstractPatchPatchDistanceFunctorT( _Scalar spatial_threshold, _Scalar angle_threshold, _Scalar scale )
        : _spatial_threshold( spatial_threshold )
        , _angle_threshold  ( angle_threshold   )
        , _scale            ( scale             )
    {}

    inline _Scalar getAngularThreshold() const { return _angle_threshold; }
    inline _Scalar getSpatialThreshold() const { return _spatial_threshold; }
    inline _Scalar getScale()            const { return _scale; }

    virtual std::string toString() const { return "AbstractPointPatchDistanceFunctorT"; }

    // relay spatial query to template class _SpatialpointPatchDistanceFunctorT
    template <class _PointT /*= typename _PointContainerT::value_type*/, class _PatchT, class _PointContainerT>
    static inline _Scalar evalSpatial( PidT const point_id, _PatchT const& patch0, _PointContainerT const& points )
    {
        _PatchT patch1; patch1.push_back( typename _PatchT::value_type( point_id, -1 ) );
        return _SpatialPatchPatchDistanceFunctorT::template eval<_PointT>( patch0, patch1, points, NULL );
    }

    // relay spatial query to template class _SpatialpointPatchDistanceFunctorT
    template <class _PointT/* = typename _PointContainerT::value_type*/, class _PatchAT, class _PatchBT, class _PointContainerT>
    static inline _Scalar evalSpatial( _PatchAT const& patch0, _PatchBT const& patch1, _PointContainerT const& points )
    {
        return _SpatialPatchPatchDistanceFunctorT::template eval<_PointT>( patch0, patch1, points, NULL );
    }

    protected:
        _Scalar _spatial_threshold, _angle_threshold, _scale;
}; //...struct AbstractPointPatchDistanceFunctorT

#if RAPTER_WITH_FULL_LINKAGE

//! \brief Distance between patches is the maximum angular distance, given the spatial distance is within threshold.
template <typename _Scalar
         , class   _SpatialPatchPatchDistanceFunctorT /*= SpatialPatchPatchMaxDistanceFunctorT<_Scalar>*/ > // Max: proper full linkage, Min: hybrid linkage (max angle, min space)
struct FullLinkagePatchPatchDistanceFunctorT : public AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>
{
        FullLinkagePatchPatchDistanceFunctorT( _Scalar spatial_threshold, _Scalar angle_threshold, _Scalar scale )
            : AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>( spatial_threshold, angle_threshold, scale ) {}

        template <class _PointT, class PatchAT, class PatchBT, class PointContainerT>
        inline _Scalar eval( PatchAT               const& patch0
                           , PatchBT               const& patch1
                           , PointContainerT       const& points
                           , _Scalar               const* current_min ) const // TODO: use current_min to stop early
        {
            typedef typename _PointT::VectorType VectorType; // concept: Eigen::Vector3f

            // get max distance from point to points in patch
            const _Scalar spatial_thresh = this->getSpatialThreshold();
            const _Scalar spatial_distance = _SpatialPatchPatchDistanceFunctorT::template eval<_PointT>( patch0, patch1, points, &spatial_thresh );

            // if close enough spatially to patch
            if ( spatial_distance > spatial_thresh )
            {
                // should not be patched together, if spatially not close enough
                return std::numeric_limits<_Scalar>::max();
            }

            // get max angular distance between line at point and lines in patch
            _Scalar max_angle( 0 );
            for ( size_t pid_id0 = 0; pid_id0 != patch0.size(); ++pid_id0 )
            {
                const int pid0  = patch0[ pid_id0 ].first;
                for ( size_t pid_id1 = 0; pid_id1 != patch1.size(); ++pid_id1 )
                {
                    const int pid1  = patch1[ pid_id1 ].first;
                    _Scalar   angle = angleInRad( points[pid0].template dir(), points[pid1].template dir() );
                    if ( angle > max_angle )
                        max_angle = angle;

                    // early exit
                    if ( current_min && (max_angle > *current_min) ) // we need the max angle to be below current_min
                        return max_angle;
                }
            }

            return max_angle;
        }

        inline _Scalar getThreshold() const { return this->getAngularThreshold(); }

        virtual std::string toString() const override { return "FullLinkagePointPatchDistanceFunctorT with " + _SpatialPatchPatchDistanceFunctorT::toString(); }
}; // ...struct FullLinkagePointPatchDistanceFunctorT

//! \brief Distance between patches is the minimum, weighted 6D squared distance.
template <typename _Scalar
         , class   _SpatialPatchPatchDistanceFunctorT /*= SpatialPatchPatchMinDistanceFunctorT<_Scalar>*/ >
struct SquaredPatchPatchDistanceFunctorT
        : public AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>
{
        SquaredPatchPatchDistanceFunctorT( _Scalar spatial_threshold, _Scalar angle_threshold, _Scalar scale )
            : AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>( spatial_threshold, angle_threshold, scale )
            , _sqr_dist_weight( (M_PI * M_PI) / (4. * scale * scale) ) {}

        template < class _PointT /*= typename PointContainerT::value_type */
                 , class PatchAT, class PatchBT
                 , class PointContainerT
                 >
        inline _Scalar eval( PatchAT                const& patch0
                           , PatchBT                const& patch1
                           , PointContainerT       const& points
                           , _Scalar               const* current_min
                           ) const
        {
            std::cout << "[" << __func__ << "]: " << "SQR EVAL" << std::endl;
            typedef typename _PointT::VectorType VectorType; // concept: Eigen::Vector3f

            // return, if not close enough spatially to patch
            const _Scalar spatial_threshold = this->getSpatialThreshold();
            const _Scalar spatial_distance  = _SpatialPatchPatchDistanceFunctorT::template eval<_PointT>( patch0, patch1, points, &spatial_threshold ); // get max distance from point to points in patch
            if ( spatial_distance > spatial_threshold )
            {
                // should not be patched together, if spatially not close enough
                return std::numeric_limits<_Scalar>::max();
            }

            // get max angular distance between line at point and lines in patch
            _Scalar max_angle( 0 );
            for ( size_t pid_id0 = 0; pid_id0 != patch0.size(); ++pid_id0 )
            {
                const int pid0 = patch0[pid_id0].first;
                for ( size_t pid_id1 = 0; pid_id1 != patch1.size(); ++pid_id1 )
                {
                    const int pid1 = patch1[pid_id1].first;
                    _Scalar   angle = angleInRad( points[pid0].template dir(), points[pid1].template dir() );
                    if ( angle > max_angle )
                        max_angle = angle;

                    // early exit
                    if ( current_min && (max_angle > *current_min) ) // we need the max angle to be below current_min
                        return max_angle;
                }
            }

            // wd = (pi^2) / (4 * scale^2);
//            static const double wd = (M_PI * M_PI) / (4. * this->getScale() * this->getScale());
//            std::cout << "sqrt(max_angle * max_angle + wd * spatial_distance * spatial_distance) = "
//                      << "sqrt(" << max_angle <<" * "<< max_angle <<" + "<< wd <<" * " << spatial_distance << " * " << spatial_distance << ") = "
//                      << "sqrt(" << max_angle * max_angle << " + " << wd << " * " << spatial_distance * spatial_distance << ") = "
//                      << sqrt( max_angle * max_angle + wd * spatial_distance * spatial_distance ) << std::endl;
            return sqrt( max_angle * max_angle + _sqr_dist_weight * spatial_distance * spatial_distance ); // 6D distance
        } // ...eval()

        inline _Scalar getThreshold() const { return this->getAngularThreshold()*this->getAngularThreshold() + _sqr_dist_weight * this->getSpatialThreshold() * this->getSpatialThreshold(); }

        virtual std::string toString() const override { return "SquaredPointPatchDistanceFunctorT with " + _SpatialPatchPatchDistanceFunctorT::toString(); }

    protected:
        _Scalar _sqr_dist_weight;

}; // ...struct SquaredPointPatchDistanceFunctorT

//! \brief Distance between patches is the angle of the representatives, given the minimum spatial distance is within threshold
template <typename _Scalar
         , class   _SpatialPatchPatchDistanceFunctorT /*= SpatialPatchPatchMaxDistanceFunctorT<_Scalar>*/ >
struct RepresentativePatchPatchDistanceFunctorT
        : public AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>
{
        RepresentativePatchPatchDistanceFunctorT( _Scalar spatial_threshold, _Scalar angle_threshold, _Scalar scale )
            : AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>( spatial_threshold, angle_threshold, scale ) {}

        template <class _PointT /*= typename PointContainerT::value_type*/, class PatchAT, class PatchBT, class PointContainerT>
        inline _Scalar eval( PatchAT               const& patch0
                           , PatchBT               const& patch1
                           , PointContainerT       const& points
                           , _Scalar               const* current_min ) const // TODO: use current_min to stop early
        {
            std::cout << "[" << __func__ << "]: " << "REPR EVAL" << std::endl;
            typedef typename _PointT::VectorType VectorType; // concept: Eigen::Vector3f

            // get max distance from point to points in patch
            const _Scalar spatial_thresh = this->getSpatialThreshold();
            const _Scalar spatial_distance = _SpatialPatchPatchDistanceFunctorT::template eval<_PointT>( patch0, patch1, points, &spatial_thresh );

            // if close enough spatially to patch
            if ( spatial_distance > spatial_thresh )
            {
                // should not be patched together, if spatially not close enough
                return std::numeric_limits<_Scalar>::max();
            }

            // get max angular distance between line at point and lines in patch
            return rapter::angleInRad( patch0.template dir(), patch1.template dir() );
        }

        inline _Scalar getThreshold() const { return this->getAngularThreshold(); }

        virtual std::string toString() const override { return "RepresentativePatchPatchDistanceFunctorT with " + _SpatialPatchPatchDistanceFunctorT::toString(); }
}; // ...struct FullLinkagePointPatchDistanceFunctorT
#endif

/*! \brief Distance is combined by spatial and angular terms according to formula: \f$ patch\_spatial\_weight^2 \cdot \frac{spat\_dist^2}{spat\_thresh} + \frac{ang\_diff^2}{ang\_thresh} < 1 \f$.
*          It's basically an ellipse in hough-space. The weight was mostly set to 0.5 \sa \ref rapter::CandidateGeneratorParams.
*   \todo  Change to actually use the direction to threshold ( ellipsoid in real space ).
*/
template < typename _Scalar
         , class    _SpatialPatchPatchDistanceFunctorT /*= SpatialPatchPatchSingleDistanceFunctorT<_Scalar> */>
struct RepresentativeSqrPatchPatchDistanceFunctorT
        : public AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>
{
        RepresentativeSqrPatchPatchDistanceFunctorT( _Scalar spatial_threshold, _Scalar angle_threshold, _Scalar scale, _Scalar dist_weight )
            : AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>( spatial_threshold, angle_threshold, scale )
            , _sqr_dist_weight   ( dist_weight       * dist_weight       )
            , _sqr_spatial_thresh( spatial_threshold * spatial_threshold )
            , _sqr_ang_thresh    ( angle_threshold   * angle_threshold   ) {}


        template <class _PointT, class _PatchAT, class _PatchBT, class _PointContainerT >
        inline _Scalar eval( _PatchAT               const& patch0
                           , _PatchBT               const& patch1
                           , _PointContainerT       const& points
                           , _Scalar                const* /*current_min*/ ) const
        {
            typedef typename _PointT::VectorType VectorType; // concept: Eigen::Vector3f

            // get max distance from point to points in patch
            const _Scalar spatial_distance = _SpatialPatchPatchDistanceFunctorT::template eval<_PointT>( patch0, patch1, points, NULL );
            const _Scalar ang_diff         = rapter::angleInRad( patch0.template dir(), patch1.template dir() );

            _Scalar val = _sqr_dist_weight * spatial_distance * spatial_distance / _sqr_spatial_thresh
                        +                    ang_diff         * ang_diff         / _sqr_ang_thresh;     // ellipse longer along line

            // get max angular distance between line at point and lines in patch
            return val;
        }

        //! \brief Returns the threshold to apply to the returned value of #eval(). (<= 1.)
        inline _Scalar getThreshold() const { return _Scalar(1.); }

        virtual std::string toString() const override { return "RepresentativePatchPatchDistanceFunctorT with " + _SpatialPatchPatchDistanceFunctorT::toString(); }

    protected:
        const _Scalar _sqr_dist_weight;
        const _Scalar _sqr_spatial_thresh;
        const _Scalar _sqr_ang_thresh;
}; // ...struct FullLinkagePointPatchDistanceFunctorT

} // ... namespace rapter

#endif // RAPTER_PATCHDISTANCEFUNCTORS_H
