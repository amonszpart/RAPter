#ifndef POINTPRIMITIVE_H
#define POINTPRIMITIVE_H

#include "primitive.h"
//#include "KMPoint.h"

#define SQR(a) ((a)*(a))

namespace am
{
    namespace kmeans
    {

        template <class TPoints,class TPoint>
        class PointPrimitive : public Primitive<TPoints,TPoint>
        {
            public:
                typedef typename TPoint::VectorType VectorType;
                typedef typename TPoint::Scalar     Scalar;

                PointPrimitive()
                    : Primitive<TPoints,TPoint>()
                    , _pnt( TPoint::VectorType::Zero() )
                {}
                PointPrimitive( VectorType pnt )
                    : Primitive<TPoints,TPoint>()
                    , _pnt( pnt )
                {}

                virtual ~PointPrimitive()
                {}

                virtual int
                recalculate( TPoints            const& points
                             , std::vector<int> const& labels
                             , int                     primitive_id );

                virtual Scalar
                getDistance( TPoint       const& point
                             , VectorType const* dim_coeffs = NULL );

                VectorType const& Pnt() const { return _pnt; }
                VectorType      & Pnt()       { return _pnt; }

            protected:
                VectorType _pnt;
        };

        template <class TPoints, class TPoint>
        int
        PointPrimitive<TPoints,TPoint>::recalculate( TPoints            const& points
                                                     , std::vector<int> const& labels
                                                     , int                     primitive_id )
        {
            // clear
            _pnt = VectorType::Zero();
            this->_size = 0;

            // sum up own labelled points
            int i = 0;
            for ( typename TPoints::const_iterator it =  points.begin();
                  /*                            */ it != points.end();
                  /*                          */ ++it, ++i )
            {
                if ( labels[i] == primitive_id )
                {
                    _pnt += it->pos();
                    ++(this->_size);
                }
            }

            // average
            if ( this->_size )
            {
                _pnt /= static_cast<Scalar>(this->_size);
            }

            return EXIT_SUCCESS;
        }

        template <class TPoints,class TPoint>
        typename TPoint::Scalar
        PointPrimitive<TPoints,TPoint>::getDistance( TPoint       const& point
                                                     , VectorType const* dim_coeffs )
        {
            VectorType diff = point.pos() - _pnt;
            if ( dim_coeffs )
                for ( int dim = 0; dim != TPoint::Dim; ++dim )
                    diff( dim ) *= dim_coeffs->operator()( dim );

            return sqrt( diff.dot(diff) );
        }
    } // nskmeans
} // nsam

#undef SQR

#endif // POINTPRIMITIVE_H
