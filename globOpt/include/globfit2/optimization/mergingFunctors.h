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

        //    std::cout << "testing (" << l0.getTag(_LineT::TAGS::GID )     << ","
        //                             << l0.getTag(_LineT::TAGS::DIR_GID ) << ")"
        //              << " vs. ("    << l1.getTag(_LineT::TAGS::GID ) << ","
        //                             << l1.getTag(_LineT::TAGS::DIR_GID ) << ")\t" << std::endl;

        // we don't merge when both group id and direction id are identical
        if (l0.getTag(_LineT::TAGS::DIR_GID ) == l1.getTag(_LineT::TAGS::DIR_GID ) &&
                l0.getTag(_LineT::TAGS::GID ) == l1.getTag(_LineT::TAGS::GID ) )
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

        //const _Scalar sqScale = scale*scale;
        //const _Scalar l0SqLengthAndScale = (l0b-l0a).squaredNorm() + scale*scale;

        // check if they have the same tag
        //bool sameTag = l0.getTag(_LineT::TAGS::DIR_GID ) == l1.getTag(_LineT::TAGS::DIR_GID );

        // check if l1 is aligned to l0
        if ( /*sameTag ||*/ // exactly aligned
             ( std::abs(l0.template normal().dot(l1a-l0a)) <= scale && // check l1a-proj(l1a,l0) <= scale
               std::abs(l0.template normal().dot(l1b-l0a)) <= scale)){  // check l1b-proj(l1b,l0) <= scale

            // check if at least one l1 endpoint is projected onto l0
            const PointT l0dir = (l0b - l0a).normalized();
            const _Scalar dl0  = (l0b - l0a).norm() + scale;
            const _Scalar dl1a = l0dir.dot(l1a-l0a);
            const _Scalar dl1b = l0dir.dot(l1b-l0a);

            if ((dl1a >= -scale && dl1a <= dl0) ||
                    (dl1b >= -scale && dl1b <= dl0))
                return true;
        }

        // check if l0 is aligned to l1
        if ( /*sameTag ||*/
             ( std::abs(l1.template normal().dot(l0a-l1a)) <= scale &&  // check l0a-proj(l0a,l0) <= scale
               std::abs(l1.template normal().dot(l0b-l1a)) <= scale )){ // check l0b-proj(l0b,l0) <= scale

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
private:
    template <class _PlaneT, class _PointContainerT, typename _Scalar>
    inline bool directionnalEval(
              _PointContainerT const& extrema0
            , _PlaneT const& p0
            , _PointContainerT const& extrema1
            , _PlaneT const& p1
            , _Scalar scale) const {
        typedef typename _PointContainerT::value_type PointT;

        // We use the following algorithm:
        // 1. Compute p0 local frame, where y is p0 normal direction.
        // 2. Express the two planes extrema in the local frame of p0.
        // 3. Check if p1 extrema y coordinate are all included in [-scale,scale]
        // 4. Check if at least one extrema is included in the finite plane p0 (according to its extrema)


        // 1. Compute p0 local frame
        //  - build a local copy p0local centered in 0 and oriented with updirection=y
        //  - compute the rotation between p0 and p0local
        //  - compute p0 gravity center
        PointT p0normal      ((extrema0[1]-extrema0[0]).normalized().cross((extrema0[2]-extrema0[1])).normalized());
        PointT p0localnormal (_Scalar(0.), _Scalar(1.), _Scalar(0.));
        //Eigen::Quaternion<_Scalar> q = Eigen::Quaternion<_Scalar>::FromTwoVectors(p0.normal(), p0local.normal());
        Eigen::Quaternion<_Scalar> q = Eigen::Quaternion<_Scalar>::FromTwoVectors(p0normal, p0localnormal);

        PointT center (PointT::Zero());
        std::for_each(extrema0.begin(), extrema0.end(), [&center] (const PointT& p){ center+=p; });
        center /= _Scalar(extrema0.size());

        // 2. Express the two planes extrema in the local frame of p0 (extremas of p0 will
        // be computed later if required
        // We use a lambda expression to first translate and then rotate the extremas
        auto transform = [&center, &q] (PointT& p) { p = q*(p-center); };
        _PointContainerT extrema1local = extrema1;
        std::for_each(extrema1local.begin(), extrema1local.end(), transform);

        // 3. Check if p1 extrema y coordinate are all included in [-scale,scale]
        auto inScaleBound = [&scale] (PointT& p){ return std::abs(p(1)) > scale; };
        typename _PointContainerT::iterator it = std::find_if(extrema1local.begin(), extrema1local.end(), inScaleBound);
        if (it != extrema1local.end()) {
            //std::cerr << "Exiting at test 1" << std::endl;
            return false;
        }
        //std::cerr << "Test 1 passed" << std::endl;

        // 4. Check if at least one extrema is included in the finite plane p0 (according to its extrema)
        _PointContainerT extrema0local = extrema0;
        std::for_each(extrema0local.begin(), extrema0local.end(), transform);


        // The computation of the local frame direction (f1, f2) is hardcoded wrt to the
        // order of the extrema computed in PlanePrimitive::getExtent
        //
        // Here we compute f1, f2, directions of the normal frame (1 coordinate shoud be = 0, we work in the plane frame)
        // and hw, hh the dimensions of the planes is those directions, extended by scale.
        // center is the middle of the finite plane
        PointT f1 = (extrema0local[1]-extrema0local[0]);
        PointT f2 = (extrema0local[2]-extrema0local[1]);
        _Scalar hw = scale + f1.norm() / _Scalar(2.); f1.normalize();
        _Scalar hh = scale + f2.norm() / _Scalar(2.); f2.normalize();


        PointT centerlocal (PointT::Zero());
        std::for_each(extrema0local.begin(), extrema0local.end(), [&centerlocal] (PointT& p){ centerlocal+=p; });
        centerlocal /= _Scalar(extrema0local.size());
        //if (centerlocal.norm() >= 10.e-6)
        //	std::cout << "problem here " << centerlocal.transpose() << std::endl;

        // here we project p over the two basis axis and check the norm of the projected vector
        // is < to the plane dimensions (hh,hw)
        auto isInFinitePlane = [&f1, &f2, &hw, &hh] (const PointT& p){
            return  (std::abs((p).dot(f1)) <= hw) &&
                    (std::abs((p).dot(f2)) <= hh);
        };
        // run the test on the coordinates of the second plane, expressed in the coordinate
        // system of the first. We merge is at least one of this coordinate is inside the
        // plane.
        it = std::find_if(extrema1local.begin(), extrema1local.end(), isInFinitePlane);
        return  it != extrema1local.end(); // true if we found someone
    }

public:
    template <class _PlaneT, class _PointContainerT, typename _Scalar>
    inline bool eval(
              _PointContainerT const& extrema0
            , _PlaneT const& p0
            , _PointContainerT const& extrema1
            , _PlaneT const& p1
            , _Scalar scale) const {
        typedef typename _PointContainerT::value_type PointT;

        // we don't merge when both group id and direction id are identical
        if (p0.getTag(_PlaneT::TAGS::DIR_GID ) == p1.getTag(_PlaneT::TAGS::DIR_GID ) &&
            p0.getTag(_PlaneT::TAGS::GID )     == p1.getTag(_PlaneT::TAGS::GID ) )
            return false;


        // compute patch with biggest area
        _Scalar a1 = (extrema0[1]-extrema0[0]).norm() * (extrema0[2]-extrema0[1]).norm();
        _Scalar a2 = (extrema1[1]-extrema1[0]).norm() * (extrema1[2]-extrema1[1]).norm();

        if (a1 > a2) return directionnalEval(extrema0, p0, extrema1, p1, scale);
        if (a1 < a2) return directionnalEval(extrema1, p1, extrema0, p0, scale);

        std::cout << "Patches with the same area: " << std::endl;
        std::cout << "(" << p0.getTag(_PlaneT::TAGS::GID )     << ","
                  << p0.getTag(_PlaneT::TAGS::DIR_GID ) << ")"
                  << " vs. ("    << p1.getTag(_PlaneT::TAGS::GID ) << ","
                  << p1.getTag(_PlaneT::TAGS::DIR_GID ) << ")\t" << std::endl;

        // same area (few chances...)
        return directionnalEval(extrema0, p0, extrema1, p1, scale) ||
               directionnalEval(extrema1, p1, extrema0, p0, scale);

        return false;
    }
};

} // namespace GF2


#endif // MERGINGFUNCTORS_H
