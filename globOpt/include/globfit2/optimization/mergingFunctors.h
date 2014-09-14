#ifndef MERGINGFUNCTORS_H
#define MERGINGFUNCTORS_H

#include <Eigen/Geometry>

namespace GF2 {
struct DecideMergeLineFunctor {
    template <class _LineT, class _PointContainerT, typename _Scalar>
    inline bool eval(
              _PointContainerT const& extrema0
            , _LineT const& l0
            , _PointContainerT const& extrema1
            , _LineT const& l1
            , _Scalar scale) const {
        typedef typename _PointContainerT::value_type PointT;

        //    std::cout << "testing (" << l0.getTag(_LineT::GID )     << ","
        //                             << l0.getTag(_LineT::DIR_GID ) << ")"
        //              << " vs. ("    << l1.getTag(_LineT::GID ) << ","
        //                             << l1.getTag(_LineT::DIR_GID ) << ")\t" << std::endl;

        // we don't merge when both group id and direction id are identical
        if (l0.getTag(_LineT::DIR_GID ) == l1.getTag(_LineT::DIR_GID ) &&
                l0.getTag(_LineT::GID ) == l1.getTag(_LineT::GID ) )
            return false;


        // We have two lines l0 and l1, respectively defined by
        // l0a-l0b and l1a-l1b
        //  A. Check if both lines have the same direction id. If not,
        //      1. We project l0a and l0b to l1 (orthogonally to l1), and compute the norm of l0a-proj(l0a,l1) and l0b-proj(l0b,l1)
        //      2. Compute the other way
        //      3. Check if the two lines are almost aligned (all norms are smaller than the scale parameters)
        //  B. If one way is correct, project second line end points to the first one
        //     and check at least one end point is projected the segment
        const PointT & l0a = extrema0[0];
        const PointT & l0b = extrema0[1];
        const PointT & l1a = extrema1[0];
        const PointT & l1b = extrema1[1];

        const PointT n0 = l0.template normal<_Scalar>();

        //const _Scalar sqScale = scale*scale;
        //const _Scalar l0SqLengthAndScale = (l0b-l0a).squaredNorm() + scale*scale;

        // check if they have the same tag
        //bool sameTag = l0.getTag(_LineT::DIR_GID ) == l1.getTag(_LineT::DIR_GID );

        // check if l1 is aligned to l0
        if ( /*sameTag ||*/ // exactly aligned
             ( std::abs(n0.dot(l1a-l0a)) <= scale && // check l1a-proj(l1a,l0) <= scale
               std::abs(n0.dot(l1b-l0a)) <= scale)){  // check l1b-proj(l1b,l0) <= scale

            // check if at least one l1 endpoint is projected onto l0
            const PointT l0dir = (l0b - l0a).normalized();
            const _Scalar dl0  = (l0b - l0a).norm() + scale;
            const _Scalar dl1a = l0dir.dot(l1a-l0a);
            const _Scalar dl1b = l0dir.dot(l1b-l0a);

            if ((dl1a >= -scale && dl1a <= dl0) ||
                    (dl1b >= -scale && dl1b <= dl0))
                return true;
        }

        const PointT n1 = l1.template normal<_Scalar>();

        // check if l0 is aligned to l1
        if ( /*sameTag ||*/
             ( std::abs(n1.dot(l0a-l1a)) <= scale &&  // check l0a-proj(l0a,l0) <= scale
               std::abs(n1.dot(l0b-l1a)) <= scale )){ // check l0b-proj(l0b,l0) <= scale

            // check if at least one l0 endpoint is projected onto l1
            const PointT l1dir = (l1b - l1a).normalized();
            const _Scalar dl1  = (l1b - l1a).norm() + scale;
            const _Scalar dl0a = l1dir.dot(l0a-l1a);
            const _Scalar dl0b = l1dir.dot(l0b-l1a);

            if ((dl0a >= -scale && dl0a <= dl1) ||
                    (dl0b >= -scale && dl0b <= dl1))
                return true;
        }

        return false;
    }
};


struct DecideMergePlaneFunctor {
    template <class _PlaneT, class _PointContainerT, typename _Scalar>
    inline bool eval(
              _PointContainerT const& /*extrema0*/
            , _PlaneT const& p0
            , _PointContainerT const& /*extrema1*/
            , _PlaneT const& p1
            , _Scalar /*scale*/) const {
        typedef typename _PointContainerT::value_type PointT;

        //    std::cout << "testing (" << p0.getTag(_LineT::GID )     << ","
        //                             << p0.getTag(_LineT::DIR_GID ) << ")"
        //              << " vs. ("    << p1.getTag(_LineT::GID ) << ","
        //                             << p1.getTag(_LineT::DIR_GID ) << ")\t" << std::endl;

        // we don't merge when both group id and direction id are identical
        if (p0.getTag(_PlaneT::DIR_GID ) == p1.getTag(_PlaneT::DIR_GID ) &&
            p0.getTag(_PlaneT::GID )     == p1.getTag(_PlaneT::GID ) )
            return false;

        // We use the following algorithm:
        // 1. Compute p0 local frame
        // 2. Express the two planes in the frame of p0, where y is the normal direction.
        // 3. Check if p1 extrema y coordinate are all included in [-scale,scale]
        // 4. Check if at least one extrema is included in the finite plane p0 (according to its extrema)

        // \todo


        // 1. Compute p0 local frame
        // We build a local copy of the plane


        return false;
    }
};

} // namespace GF2


#endif // MERGINGFUNCTORS_H
