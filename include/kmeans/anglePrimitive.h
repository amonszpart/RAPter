#ifndef ANGLEPRIMITIVE_H
#define ANGLEPRIMITIVE_H

#include "pointPrimitive.h"

#define SQR(a) ((a)*(a))

namespace am
{
    namespace kmeans
    {

        template <class TPoints,class TPoint>
        class AnglePrimitive : public PointPrimitive<TPoints,TPoint>
        {
            public:
                typedef typename TPoint::VectorType VectorType;
                typedef typename TPoint::Scalar     Scalar;

                AnglePrimitive(  std::vector<int> const& angle_dims )
                    : PointPrimitive<TPoints,TPoint>()
                    , _angle_dims( angle_dims )
                {}

                AnglePrimitive( VectorType pnt, std::vector<int> const& angle_dims )
                    : PointPrimitive<TPoints,TPoint>( pnt )
                    , _angle_dims( angle_dims )
                {}

                virtual Scalar
                getDistance( TPoint       const& point
                             , VectorType const* dim_coeffs = NULL );
            protected:
                std::vector<int> _angle_dims;

                // hide these
                AnglePrimitive( VectorType pnt ) {};
                AnglePrimitive() {};
        };

        template <class TPoints,class TPoint>
        typename TPoint::Scalar
        AnglePrimitive<TPoints,TPoint>::getDistance( TPoint       const& point
                                                     , VectorType const* dim_coeffs )
        {
            //        if ( dim_coeffs )
            //            std::cout << __func__ << "]  dim_coeffs: " << dim_coeffs->transpose()
            //                      << " angle_dims.size: " << _angle_dims.size() << "; angle_dims[0]: " << _angle_dims[0] << std::endl;
            //        else
            //            std::cout << __func__ << "] no dim_coeffs" << ", " << " angle_dims.size: " << _angle_dims.size() << "; angle_dims[0]: " << _angle_dims[0] << std::endl;
            VectorType full_diff = point.pos() - this->_pnt;
            VectorType diff( VectorType::Zero() );
            for ( int dim_id = 0; dim_id != _angle_dims.size(); ++dim_id )
            {
                const int dim = _angle_dims[ dim_id ];
                diff( dim )   = full_diff( dim );

                int i = (int)(diff(dim) / M_PI);
                //            if ( i )
                //            {
                //                std::cout << "dim_id:\t" << dim_id << std::endl;
                //                std::cout << "diff(" << dim << "):\t" << diff(dim) << std::endl;
                //                std::cout << "(int)(diff(dim) / M_PI):\t" << i << std::endl;
                //                std::cout << "diff between:\t" << point.pos()(dim) * 180./M_PI << " and " << this->_pnt(dim) * 180./M_PI << std::endl;
                //                std::cout << "diff(dim)-i*M_PI:\t" << (diff(dim) - i * M_PI) * 180./M_PI << std::endl << std::endl;
                //            }
                diff(dim) -= i * M_PI;

                if ( dim_coeffs )
                    diff( dim ) *= dim_coeffs->operator()( dim );
            }

            return sqrt( diff.dot(diff) );
        }

    } // ns kmeans

} // nsam

#undef SQR

#endif // ANGLEPRIMITIVE_H

//                std::cout << "dim_id:\t" << dim_id << std::endl;
//                std::cout << "diff(" << dim << "):\t" << diff(dim) << std::endl;
//                std::cout << "(int)(diff(dim) / M_PI_2):\t" << i << std::endl;
//                std::cout << "diff between:\t" << point.pos()(dim) * 180./M_PI << " and " << this->_pnt(dim) * 180./M_PI << std::endl;
//                std::cout << "diff(dim)-i*M_PI_2:\t" << (diff(dim) - i * M_PI_2) * 180./M_PI << std::endl << std::endl;
