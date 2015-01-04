#ifndef __GF2_ENERGYFUNCTORS_H__
#define __GF2_ENERGYFUNCTORS_H__

#include <vector>
#include "globfit2/my_types.h" // angleInRad
#include "globfit2/parameters.h"

#include <Eigen/StdVector>

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
                // Changed by Aron on 18/09/2014
                return std::abs( primitive.getDistance(point) );
                // this only works for lines:
                // return (primitive.pos() - point.pos()).cross(primitive.dir()).norm();
            }
    };

    struct SharedVolumeForPlanesWithScaleFunctor
    {
        // return 1 for perfect match, 0 for invalid
        template < class _PrimitiveT
                 , class _PointContainerT
                 , typename _Scalar>
        inline typename _PrimitiveT::Scalar eval(
                  _PointContainerT const& extrema0
                , _PrimitiveT const& p0
                , _PointContainerT const& extrema1
                , _PrimitiveT const& p1
                , _Scalar scale) const {
            typedef typename _PointContainerT::value_type PointT;

            // Compute the volume of the boxes defined by the extrema and the scale (both up and bottom)
            _Scalar v0 = (extrema0[1]-extrema0[0]).norm() * (extrema0[2]-extrema0[1]).norm() * _Scalar(2.) * scale;
            _Scalar v1 = (extrema1[1]-extrema1[0]).norm() * (extrema1[2]-extrema1[1]).norm() * _Scalar(2.) * scale;

            _PointContainerT points = extrema0;
            points.insert(points.end(), extrema1.begin(), extrema1.end());

            PointT center (PointT::Zero());
            std::for_each(points.begin(), points.end(), [&center] (const PointT& p){ center+=p; });
            center /= _Scalar(extrema0.size());

            Eigen::Matrix<_Scalar, 3, 3> cov = Eigen::Matrix<_Scalar, 3, 3>::Zero();
            std::for_each(points.begin(), points.end(), [&center, &cov] (const PointT& p){ cov+=(p-center)*(p-center).transpose(); });

            Eigen::SelfAdjointEigenSolver< Eigen::Matrix<_Scalar, 3, 3> > es;
            es.compute( cov );

            // Build a matrix n*3, where n is the number of points,
            // By multiplying with the eigenvector matrix, we get the projected distances
            // of all the points in each directions.
            // we can then extract the min/max of these distances, and compute the volume
            // of the box.
            Eigen::Matrix<_Scalar, Eigen::Dynamic, 3> relativePoints (points.size(), 3);
            for(unsigned int i = 0; i != points.size(); ++i)
                relativePoints.row(i) = (points(i) - center).transpose();

            Eigen::Matrix<_Scalar, 3, 3> projMatrix = relativePoints * es.eigenvectors();
            _Scalar volume  = (projMatrix.colwise().maxCoeff().array() - projMatrix.colwise().minCoeff().array()).prod();

            // we know that volume is necessary bigger than v0 and v1
            return _Scalar(1.) - (volume - (v0+v1)) / volume;
        }
    };

    struct SharedAreaForLinesWithScaleFunctor
    {
        // return 1 for perfect match, 0 for invalid
        template < class _PrimitiveT
                 , class _PointContainerT
                 , typename _Scalar>
        inline typename _PrimitiveT::Scalar eval(
                  _PointContainerT const& extrema0
                , _PrimitiveT const& p0
                , _PointContainerT const& extrema1
                , _PrimitiveT const& p1
                , _Scalar scale) const {
            typedef typename _PointContainerT::value_type PointT;
            typedef Eigen::Matrix<_Scalar, 2, 1> Point2DT;

            // Compute the area of the boxes defined by the extrema and the scale (above and below the curve)
            _Scalar a0 = (extrema0[1]-extrema0[0]).norm() * _Scalar(2.) * scale;
            _Scalar a1 = (extrema1[1]-extrema1[0]).norm() * _Scalar(2.) * scale;


            std::vector< Point2DT > points;
            std::for_each(extrema0.begin(), extrema0.end(), [&points] (const PointT &p){ points.push_back(p.template head<2>()); });
            std::for_each(extrema1.begin(), extrema1.end(), [&points] (const PointT &p){ points.push_back(p.template head<2>()); });

            Point2DT center (Point2DT::Zero());
            std::for_each(points.begin(), points.end(), [&center] (const Point2DT& p){ center+=p; });
            center /= _Scalar(points.size());

            Eigen::Matrix<_Scalar, 2, 2> cov = Eigen::Matrix<_Scalar, 2, 2>::Zero();
            std::for_each(points.begin(), points.end(), [&center, &cov] (const Point2DT& p){ cov+=(p-center)*(p-center).transpose(); });

            Eigen::SelfAdjointEigenSolver< Eigen::Matrix<_Scalar, 2, 2> > es;
            es.compute( cov );

            // Build a matrix n*2, where n is the number of points,
            // By multiplying with the eigenvector matrix, we get the projected distances
            // of all the points in each directions.
            // we can then extract the min/max of these distances, and compute the volume
            // of the box.
            Eigen::Matrix<_Scalar, Eigen::Dynamic, 2> relativePoints (points.size(), 2);
            for(unsigned int i = 0; i != points.size(); ++i)
                relativePoints.row(i) << points.at(i).transpose();

            Eigen::Matrix<_Scalar, Eigen::Dynamic, 2> projMatrix = relativePoints * es.eigenvectors();
            _Scalar area  = (projMatrix.colwise().maxCoeff().array() - projMatrix.colwise().minCoeff().array()).prod();

            // we know that volume is necessary bigger than v0 and v1
            return _Scalar(1.) - (area - (a0+a1)) / area;
        }
    };

    //! \brief Compute the distance to a finite line
    struct MyPointFiniteLineDistanceFunctor {
        template <class _LineT, class _PointContainerT, class _Vec3Derived>
        static inline typename _LineT::Scalar eval(
                  _PointContainerT const& extrema0
                , _LineT const& l0
                , _Vec3Derived const& q
                ) /*const*/ {
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
        static inline typename _LineT::Scalar eval(
                  _PointContainerT const& extrema0
                , _LineT const& l0
                , _Vec3Derived const& q
                ) /*const*/ {
            typedef typename _PointContainerT::value_type PointT;
            typedef typename _LineT::Scalar Scalar;

            PointT center (PointT::Zero());
            std::for_each(extrema0.begin(), extrema0.end(), [&center] (const PointT& p){ center+=p; });
            center /= Scalar(extrema0.size());

            // First we build the local frame of the finite plane using the extrema.
            // We use here a trick based on the order of these extrema (see PlanePrimitive::getExtent),
            // also used in DecideMergePlaneFunctor::eval()
            // The local frame is stored in a matrix, one column for each axis

            Eigen::Matrix<Scalar, 3,3> frame;
            frame.col(0) = (extrema0[1]-extrema0[0]).normalized();
            frame.col(1) = (extrema0[2]-extrema0[1]).normalized();
            frame.col(2) = frame.col(0).cross(frame.col(1)).normalized();

            Eigen::Matrix<Scalar, 1, 3> hsize ( // row vector
                        (extrema0[1]-extrema0[0]).norm() / Scalar(2.),
                        (extrema0[2]-extrema0[1]).norm() / Scalar(2.),
                        Scalar(0.) // a plane is flat, no thickness
                        );

            Eigen::Matrix<Scalar, 1, 3> lq = q - center; // column vector

            return (( lq * frame )                                   // project each component of lq onto the frame and get the distance (dot product)
                    .array().abs().eval().array()                    // get unsigned distance
                    <= hsize.array())                                // check if we are inside the box for each coordinate
                    .select (                                        // if test (dimension wise):
                        Eigen::Matrix<Scalar, 1, 3>::Zero(),          //   - we are inside the box, so distance is 0
                        (( lq * frame ).array().abs() - hsize.array()).eval() //   - we are outside the box, so return unsigned distance (distance to center - box size)
                        ).norm();




            return std::numeric_limits<typename _PointContainerT::value_type::Scalar>::max();
        }
    };

    template <class _PrimitiveT/*, class PointToPrimFunctor*/>
    struct MyFinitePrimitiveToFinitePrimitiveCompatFunctor
    {
        template <class _PointContainerT>
        static inline typename _PrimitiveT::Scalar eval(
                  _PointContainerT  const& extrema0
                , _PrimitiveT       const& p0
                , _PointContainerT  const& extrema1
                , _PrimitiveT       const& p1 ) /*const*/
        {
            typedef typename _PointContainerT::value_type PointT;
            typedef typename _PrimitiveT::Scalar Scalar;

            //PointToPrimFunctor functor;

            Scalar min = std::numeric_limits<Scalar>::max();
            for ( typename _PointContainerT::const_iterator it  = extrema0.begin();
                                                            it != extrema0.end  (); ++it )
            {
                //const Scalar m = functor.eval(extrema1, p1, *it);
                const Scalar m = p1.getFiniteDistance( extrema1, *it );
                if ( m < min ) min = m;
            }

            for (typename _PointContainerT::const_iterator it = extrema1.begin();
                 it != extrema1.end(); ++it){
                //const Scalar m = functor.eval(extrema0, p0, *it);
                const Scalar m = p0.getFiniteDistance( extrema0, *it );
                if(m < min) min = m;
            }

            return min;
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
        eval( const PrimitiveT &/*p1*/, const PrimitiveT &/*p2*/ )
        {
            std::cerr << "[" << __func__ << "]: " << "Abstract function, use specialization!" << std::endl;
            return std::numeric_limits<Scalar>::max();
        }

        inline std::vector<Scalar> getAngles() const { return _angles; }
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
        eval( const PrimitiveT &p1, const PrimitiveT &p2 )
        {
            Scalar diff  = MyPrimitivePrimitiveAngleFunctor::eval( p1, p2, this->_angles );

            // truncated at half degree difference
            Scalar score = std::min( Scalar(0.09341652027), Scalar(sqrt(diff)) );

            return score;
        }
    }; //...SqrtPrimitivePrimitiveEnergyFunctor

#if 1
    /*! \brief Used to setup spatially smooth pairwise cost in \ref GF2::problemSetup::formulate().
     * \tparam _PrimitiveCompFunctor Concept: \ref GF2::MyFinitePrimitiveToFinitePrimitiveCompatFunctor
     */
    template <class _FiniteFiniteDistanceFunctor, class _PointContainerT, typename _Scalar, class PrimitiveT>
    struct SpatialSqrtPrimitivePrimitiveEnergyFunctor : public AbstractPrimitivePrimitiveEnergyFunctor<_Scalar,PrimitiveT>
    {
        typedef std::vector<_Scalar> AnglesT;
        typedef std::vector< Eigen::Matrix<_Scalar,3,1> > ExtremaT;

        SpatialSqrtPrimitivePrimitiveEnergyFunctor( AnglesT              const& angles
                                                  , _PointContainerT     const& points
                                                  , _Scalar              const  scale )
            : AbstractPrimitivePrimitiveEnergyFunctor<_Scalar,PrimitiveT>( angles )
            , _points( points )
            , _scale ( scale )
            , _verbose( false )
            , _dirIdBias( 0. )
            , _truncAngle( 0. )
            , _useAngleGen( 0 )
            , _spatialWeightCoeff( 0. )
        {}

        virtual ~SpatialSqrtPrimitivePrimitiveEnergyFunctor() {};

        virtual inline _Scalar
        evalSpatial( PrimitiveT const& p1, ExtremaT const& ex1, PrimitiveT const& p2, ExtremaT const& ex2 )
        {
            // distance between finite primitives
            //_Scalar dist = _primPrimCompFunctor.template eval( ex1, p1, ex2, p2 ); // changed by Aron on 21 12 2014
            _Scalar dist = _FiniteFiniteDistanceFunctor::eval( ex1, p1, ex2, p2 );
            // (scale - distance) truncated at scale-scale = 0, normalized by scale to 0..1
            _Scalar spat_w = std::max( ProblemSetupParams<_Scalar>::spatial_weight_distance * _scale - dist
                                     , _Scalar(0.) )
                                     / _scale;

            return _spatialWeightCoeff * spat_w;
        }

        virtual inline _Scalar
        eval( PrimitiveT const& p1, ExtremaT const& ex1, PrimitiveT const& p2, ExtremaT const& ex2
            , AnglesT const& allowedAngles, _Scalar *idealAngle = NULL, _Scalar *spatialWeightOut = NULL  )
        {
            _Scalar score        ( 0. );
            _Scalar spatialWeight( 0. );

            if ( p1.getTag(PrimitiveT::TAGS::DIR_GID) != p2.getTag(PrimitiveT::TAGS::DIR_GID)          )
            {
                if (     !_useAngleGen
                     || (p1.getTag(PrimitiveT::TAGS::GEN_ANGLE) == p2.getTag(PrimitiveT::TAGS::GEN_ANGLE)) )
                {
                    //_Scalar diff  = MyPrimitivePrimitiveAngleFunctor::eval( p1, p2, this->_angles );
                    int     idealAngleId = 0;
                    _Scalar diff         = MyPrimitivePrimitiveAngleFunctor::eval( p1, p2, allowedAngles, &idealAngleId );
                    if ( idealAngle )
                        *idealAngle = allowedAngles[ idealAngleId ];

                    if (_truncAngle != _Scalar(0.) )
                        score = angleFunction( std::min(diff, _truncAngle) );
                }
                else
                {
                    if ( idealAngle )
                        *idealAngle = _Scalar( 0. );
                    //if ( spatialWeightOut )
                    //    spatialWeight = evalSpatial( p1, ex1, p2, ex2 );
                }

                spatialWeight = evalSpatial( p1, ex1, p2, ex2 );
                if ( spatialWeight > _Scalar(0.) )
                    score += _spatialWeightCoeff;
            }

            score += _dirIdBias;

            // output
            if ( spatialWeightOut )
                *spatialWeightOut = spatialWeight;

            return score;
        } //...eval()

            inline void setDirIdBias( _Scalar dirIdBias ) { _dirIdBias = dirIdBias; }
            inline void setTruncAngle( _Scalar truncAngle ) { _truncAngle = truncAngle; }
            inline void setUseAngleGen( int useAngleGen ) { _useAngleGen = useAngleGen; }
            inline int getUseAngleGen() { return _useAngleGen; }
            inline void setSpatialWeightCoeff( _Scalar spatialWeightCoeff ) { _spatialWeightCoeff = spatialWeightCoeff; }
            inline _Scalar angleFunction( _Scalar angleDiff ) { return sqrt(angleDiff); }

            bool                           _verbose;
        protected:
            _FiniteFiniteDistanceFunctor   _primPrimCompFunctor;
            _PointContainerT        const& _points;
            _Scalar                        _scale;
            _Scalar                        _dirIdBias;
            _Scalar                        _truncAngle;
            int                            _useAngleGen;
            _Scalar                        _spatialWeightCoeff;
    }; //...SpatialSqrtPrimitivePrimitiveEnergyFunctor
#endif

    template <typename Scalar, class PrimitiveT>
    struct CExpPrimitivePrimitiveEnergyFunctor : public AbstractPrimitivePrimitiveEnergyFunctor<Scalar,PrimitiveT>
    {
#if __cplusplus > 199711L
        using AbstractPrimitivePrimitiveEnergyFunctor<Scalar,PrimitiveT>::AbstractPrimitivePrimitiveEnergyFunctor;
#endif
        virtual ~CExpPrimitivePrimitiveEnergyFunctor() {};

        virtual inline Scalar
        eval( const PrimitiveT &p1, const PrimitiveT &p2 )
        {
            Scalar diff  = MyPrimitivePrimitiveAngleFunctor::eval( p1, p2, this->_angles );
            Scalar score = diff * diff * diff;
            score *= score;

            return score;
        }
    }; //...CExpPrimitivePrimitiveEnergyFunctor

} // ... ns GF2


#endif // __GF2_ENERGYFUNCTORS_H__
