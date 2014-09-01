#ifndef GF2_PATCHDISTANCEFUNCTORS_H
#define GF2_PATCHDISTANCEFUNCTORS_H

#include "globfit2/my_types.h" // angleInRad

namespace GF2
{

//! \brief SpatialPatchPatchMaxDistanceFunctorT     Spatial patch-patch distance is the maximum distance between any two members (farthest two points).
template <typename _Scalar>
struct SpatialPatchPatchMaxDistanceFunctorT
{
        // spatial distanceToPatch
        // Concept: PointContainerT = std::vector<PointPrimitive>, PatchT = std::vector<int>, PointT = PointPrimitive
        template <class _PointContainerT, class _PatchT, class _PointT = typename _PointContainerT::value_type>
        static inline _Scalar eval( _PatchT const& patch0, _PatchT const& patch1, _PointContainerT const& points, _Scalar const* cut_off )
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

//! \brief SpatialPatchPatchMinDistanceFunctorT     Spatial patch-patch distance is the minimum distance between any two members (closest two points).
template <typename _Scalar>
struct SpatialPatchPatchMinDistanceFunctorT
{
        // Concept: PointContainerT = std::vector<PointPrimitive>, PatchT = std::vector<int>, PointT = PointPrimitive
        template <class _PointContainerT, class _PatchT, class _PointT = typename _PointContainerT::value_type>
        static inline _Scalar eval( _PatchT const& patch0, _PatchT const& patch1, _PointContainerT const& points, _Scalar const* /*cut_off*/ )
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

//_________________________________________________________________________________________________________________________________________________________________________________

//! \brief AbstractPatchPatchDistanceFunctorT       Calculates and merges distance in spatial and angular domain between two patches.
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
    template <class _PointContainerT, class _PatchT, class _PointT = typename _PointContainerT::value_type>
    static inline _Scalar evalSpatial( int const point_id, _PatchT const& patch0, _PointContainerT const& points )
    {
        _PatchT patch1 = { typename _PatchT::value_type( point_id, -1 ) };
        return _SpatialPatchPatchDistanceFunctorT::eval( patch0, patch1, points, NULL );
    }

    // relay spatial query to template class _SpatialpointPatchDistanceFunctorT
    template <class _PointContainerT, class _PatchT, class _PointT = typename _PointContainerT::value_type>
    static inline _Scalar evalSpatial( _PatchT const& patch0, _PatchT const& patch1, _PointContainerT const& points )
    {
        return _SpatialPatchPatchDistanceFunctorT::eval( patch0, patch1, points, NULL );
    }

    protected:
    _Scalar _spatial_threshold, _angle_threshold, _scale;
}; //...struct AbstractPointPatchDistanceFunctorT

//! \brief FullLinkagePatchPatchDistanceFunctorT    Distance between patches is the maximum angular distance, given the spatial distance is within threshold.
template <typename _Scalar
         , class   _SpatialPatchPatchDistanceFunctorT = SpatialPatchPatchMaxDistanceFunctorT<_Scalar> > // Max: proper full linkage, Min: hybrid linkage (max angle, min space)
struct FullLinkagePatchPatchDistanceFunctorT : public AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>
{
        using AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>::AbstractPatchPatchDistanceFunctorT; //!< \brief Constructor relay

        template <class PatchT, class PointContainerT, class PointT = typename PointContainerT::value_type >
        inline _Scalar eval( PatchT                const& patch0
                           , PatchT                const& patch1
                           , PointContainerT       const& points
                           , _Scalar               const* current_min ) const // TODO: use current_min to stop early
        {
            typedef typename PointT::VectorType VectorType; // concept: Eigen::Vector3f

            // get max distance from point to points in patch
            const _Scalar spatial_thresh = this->getSpatialThreshold();
            const _Scalar spatial_distance = _SpatialPatchPatchDistanceFunctorT::eval( patch0, patch1, points, &spatial_thresh );

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

        virtual std::string toString() const override { return "FullLinkagePointPatchDistanceFunctorT with " + _SpatialPatchPatchDistanceFunctorT::toString(); }
}; // ...struct FullLinkagePointPatchDistanceFunctorT

//! \brief SquaredPatchPatchDistanceFunctorT        Distance between patches is the minimum, weighted 6D squared distance.
template <typename _Scalar
         , class   _SpatialPatchPatchDistanceFunctorT = SpatialPatchPatchMinDistanceFunctorT<_Scalar> >
struct SquaredPatchPatchDistanceFunctorT : public AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>
{
        //! \brief Parent relay
        using AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>::AbstractPatchPatchDistanceFunctorT;

        template <class PatchT
                 , class PointContainerT
                 , class PointT = typename PointContainerT::value_type >
        inline _Scalar eval( PatchT                const& patch0
                           , PatchT                const& patch1
                           , PointContainerT       const& points
                           , _Scalar               const* current_min
                           ) const
        {
            typedef typename PointT::VectorType VectorType; // concept: Eigen::Vector3f

            // return, if not close enough spatially to patch
            const _Scalar spatial_threshold = this->getSpatialThreshold();
            const _Scalar spatial_distance = _SpatialPatchPatchDistanceFunctorT::eval( patch0, patch1, points, &spatial_threshold ); // get max distance from point to points in patch
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
            static const double wd = (M_PI * M_PI) / (4. * this->getScale() * this->getScale());
//            std::cout << "sqrt(max_angle * max_angle + wd * spatial_distance * spatial_distance) = "
//                      << "sqrt(" << max_angle <<" * "<< max_angle <<" + "<< wd <<" * " << spatial_distance << " * " << spatial_distance << ") = "
//                      << "sqrt(" << max_angle * max_angle << " + " << wd << " * " << spatial_distance * spatial_distance << ") = "
//                      << sqrt( max_angle * max_angle + wd * spatial_distance * spatial_distance ) << std::endl;
            return sqrt( max_angle * max_angle + wd * spatial_distance * spatial_distance ); // 6D distance
        } // ...eval()

        virtual std::string toString() const override { return "SquaredPointPatchDistanceFunctorT with " + _SpatialPatchPatchDistanceFunctorT::toString(); }

}; // ...struct SquaredPointPatchDistanceFunctorT

//! \brief RepresentativePatchPatchDistanceFunctorT Distance between patches is the angle of the representatives, given the minimum spatial distance is within threshold
template <typename _Scalar
         , class   _SpatialPatchPatchDistanceFunctorT = SpatialPatchPatchMaxDistanceFunctorT<_Scalar> >
struct RepresentativePatchPatchDistanceFunctorT : public AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>
{
        using AbstractPatchPatchDistanceFunctorT<_Scalar, _SpatialPatchPatchDistanceFunctorT>::AbstractPatchPatchDistanceFunctorT; //!< \brief Constructor relay

        template <class PatchT, class PointContainerT, class PointT = typename PointContainerT::value_type >
        inline _Scalar eval( PatchT                const& patch0
                           , PatchT                const& patch1
                           , PointContainerT       const& points
                           , _Scalar               const* current_min ) const // TODO: use current_min to stop early
        {
            typedef typename PointT::VectorType VectorType; // concept: Eigen::Vector3f

            // get max distance from point to points in patch
            const _Scalar spatial_thresh = this->getSpatialThreshold();
            const _Scalar spatial_distance = _SpatialPatchPatchDistanceFunctorT::eval( patch0, patch1, points, &spatial_thresh );

            // if close enough spatially to patch
            if ( spatial_distance > spatial_thresh )
            {
                // should not be patched together, if spatially not close enough
                return std::numeric_limits<_Scalar>::max();
            }

            // get max angular distance between line at point and lines in patch
            return GF2::angleInRad( patch0.template dir(), patch1.template dir() );
        }

        virtual std::string toString() const override { return "RepresentativePatchPatchDistanceFunctorT with " + _SpatialPatchPatchDistanceFunctorT::toString(); }
}; // ...struct FullLinkagePointPatchDistanceFunctorT

} // ... namespace GF2
#endif // GF2_PATCHDISTANCEFUNCTORS_H
