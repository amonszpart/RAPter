#ifndef __GF2_ENERGYFUNCTORS_H__
#define __GF2_ENERGYFUNCTORS_H__

#include <vector>
#include "globfit2/my_types.h" // angleInRad

namespace GF2
{
    //! \brief Functor with eval function; Takes a point and a primitive, and returns their distance.
    //! \todo Move this to the primitive class, since this is dependent on the primitive (i.e. a cylinder won't generalize to this concept)
    struct MyPointPrimitiveDistanceFunctor
    {
            //! \brief                  \copydoc MyPointPrimitiveDistanceFunctor
            //! \tparam Scalar          Scalar type to calculate angle in. Concept: float.
            //! \tparam PointT          Point wrapper type with pos() function implemented. Concept: PointPrimitive.
            //! \tparam PrimitiveT      A wrapper type for a primitive with pos() and dir() functions implemented. Concepts: LinePrimitive, LinePrimitive2, PlanePrimitive.
            //! \param[in] point        Subject of the distance to primitive calculation.
            //! \param[in] primitive    The primitive, from which the distance is to be calculated.
            //! \return                 The distance between the point and the primitive.
            template <typename Scalar, class PointT, class PrimitiveT>
            static inline Scalar
            eval( PointT const& point, PrimitiveT const& primitive )
            {
                return (primitive.pos() - point.pos()).cross(primitive.dir()).norm();
            }
    };

    //! \brief Compute the distance to a finite line
    struct MyPointFiniteLineDistanceFunctor {
        template <class _LineT, class _PointContainerT, class _Vec3Derived>
        inline typename _LineT::Scalar eval(
                  _PointContainerT const& extrema0
                , _LineT const& l0
                , _Vec3Derived const& q
                ) const {
            typedef typename _PointContainerT::value_type PointT;
            typedef typename _LineT::Scalar Scalar;


            // The algorithm is very similar to the one used in DecideMergeLineFunctor::eval,
            // except that we need to compute the distance between the line l0 and q (and not to l1 extents)
            // There are two cases:
            //   - q is projected othogonaly to the finite line l0, in that case we return the orthogonal direction
            //   - q is projected out of the line, so we return the minimum distance to the two extents

            // check if q is in the band aligned with l0, thickness=scale
            const PointT & l0a = extrema0[0];
            const PointT & l0b = extrema0[1];

            // Get line direction using extents
            const PointT l0dir = (l0b - l0a).normalized();
            // Compute dq, the parametric position of q on the line
            const Scalar dq   = l0dir.dot(q-l0a);

            // if orthogonally projected onto the finite line
            if ((dq >= Scalar(0.) && dq <= (l0b - l0a).norm())){
                return std::abs(l0.template normal().dot(q-l0.pos())); // orthogonal direction to the line
            }
            return std::min( (q-l0a).norm(),(q-l0b).norm()); // minimum distance the the line extremas
        }
    };

    //! \brief Compute the distance to a finite line
    struct MyPointFinitePlaneDistanceFunctor {
        template <class _LineT, class _PointContainerT, class _Vec3Derived>
        inline bool eval(
                  _PointContainerT const& extrema0
                , _LineT const& l0
                , _Vec3Derived const& q
                ) const {
            return std::numeric_limits<typename _PointContainerT::value_type::Scalar>::max();
        }
    };

    //! \brief Functor with eval function; Takes two primitives, calculates their angle,
    //!        and returns the abs difference to the closest angle provided in the angles parameter.
    struct MyPrimitivePrimitiveAngleFunctor
    {

        //! \brief                      \copydoc MyPrimitivePrimitiveAngleFunctor
        //! \tparam Scalar              Scalar type to calculate angle in. Concept: float.
        //! \tparam PrimitiveT          A wrapper type for a primitive with dir() function implemented. Concepts: LinePrimitive, LinePrimitive2, PlanePrimitive.
        //! \param[in] p1               First primitive for angle calculation.
        //! \param[in] p2               Second primitive for angle calculation.
        //! \param[in] angles           Vector of angles to take absolute differences from. This function looks for the smallest distance to one of these.
        //! \param[out] closest_angle   Pointer to output the angle from \p angles that was the closest to the return value.
        //! \return                     The absolute angle difference to the closest angle in \p angles.
        template <typename Scalar, class PrimitiveT>
        static inline Scalar
        eval( PrimitiveT            const& p1
            , PrimitiveT            const& p2
            , std::vector<Scalar>   const& angles
            , int                        * closest_angle_id = NULL )
        {
            // angle
            Scalar angle = GF2::angleInRad( p1.dir(), p2.dir() );
            // check nan
            if ( angle != angle )   angle =  Scalar(0);
            // normalize to 0..180
            while ( angle > M_PI )  angle -= M_PI;

            // track closest angle
            Scalar min_angle = std::numeric_limits<Scalar>::max();
            // calculate minimum distance from all input angles
            for ( size_t i = 0; i != angles.size(); ++i )
            {
                // calc distance
                Scalar diff = std::abs( angles[i] - angle );

                // select min
                if ( diff < min_angle )
                {
                    min_angle = diff;
                    // output, if requested
//                    if ( closest_angle )
//                        *closest_angle = angles[i];
                    if ( closest_angle_id )
                        *closest_angle_id = i;
                } //...select min
            } //...for all angles

            return min_angle;
        } //...eval()
    }; //...MyPrimitivePrimitiveAngleFunctor

    template <typename Scalar, class PrimitiveT>
    struct AbstractPrimitivePrimitiveEnergyFunctor
    {
        AbstractPrimitivePrimitiveEnergyFunctor( std::vector<Scalar> angles )
            : _angles( angles ) {}
        virtual ~AbstractPrimitivePrimitiveEnergyFunctor() {}

        virtual inline Scalar
        eval( PrimitiveT /*p1*/, PrimitiveT /*p2*/ )
        {
            std::cerr << "[" << __func__ << "]: " << "Abstract function, use specialization!" << std::endl;
            return std::numeric_limits<Scalar>::max();
        }

        protected:
            std::vector<Scalar> _angles;
    }; //...AbstractPrimitivePrimitiveEnergyFunctor

    template <typename Scalar, class PrimitiveT>
    struct SqrtPrimitivePrimitiveEnergyFunctor : public AbstractPrimitivePrimitiveEnergyFunctor<Scalar,PrimitiveT>
    {
#if __cplusplus > 199711L
        using AbstractPrimitivePrimitiveEnergyFunctor<Scalar,PrimitiveT>::AbstractPrimitivePrimitiveEnergyFunctor;
#endif
        virtual ~SqrtPrimitivePrimitiveEnergyFunctor() {};

        virtual inline Scalar
        eval( PrimitiveT p1, PrimitiveT p2 )
        {
            Scalar diff  = MyPrimitivePrimitiveAngleFunctor::eval( p1, p2, this->_angles );
            Scalar score = sqrt( diff );

            return score;
        }
    }; //...SqrtPrimitivePrimitiveEnergyFunctor

    template <typename Scalar, class PrimitiveT>
    struct CExpPrimitivePrimitiveEnergyFunctor : public AbstractPrimitivePrimitiveEnergyFunctor<Scalar,PrimitiveT>
    {
#if __cplusplus > 199711L
        using AbstractPrimitivePrimitiveEnergyFunctor<Scalar,PrimitiveT>::AbstractPrimitivePrimitiveEnergyFunctor;
#endif
        virtual ~CExpPrimitivePrimitiveEnergyFunctor() {};

        virtual inline Scalar
        eval( PrimitiveT p1, PrimitiveT p2 )
        {
            Scalar diff  = MyPrimitivePrimitiveAngleFunctor::eval( p1, p2, this->_angles );
            Scalar score = diff * diff * diff;
            score *= score;

            return score;
        }
    }; //...CExpPrimitivePrimitiveEnergyFunctor

} // ... ns GF2


#endif // __GF2_ENERGYFUNCTORS_H__
