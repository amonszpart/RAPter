#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <vector>
#include "KMPoint.h"

namespace am
{
    namespace kmeans
    {
        // pre-declare
        template <class TPoints, class TPoint> class PointPrimitive;

        template <class TPoints,class TPoint>
        class Primitive
        {
            public:
                typedef std::vector<Primitive<TPoints,TPoint>*> Vec;

                virtual int
                recalculate( TPoints            const& points
                             , std::vector<int> const& labels
                             , int                     primitive_id ) = 0;

                virtual typename TPoint::Scalar
                getDistance( TPoint                        const& point
                             , typename TPoint::VectorType const* dim_coeffs = NULL ) = 0;

                virtual ~Primitive() {}
                PointPrimitive<TPoints,TPoint>* asPointPrimitive() { return reinterpret_cast<PointPrimitive<TPoints,TPoint>*>( this ); }
                int size() { return _size; }

            protected:
                int _size;
        };

    } // ns kmeans
} // ns am

#endif // PRIMITIVE_H
