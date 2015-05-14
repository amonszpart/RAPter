#ifndef HOUGHLINE_H
#define HOUGHLINE_H

#include <fstream>

#include <Eigen/Dense>
#include "primitives/primitive.h"

#ifdef GF2_USE_PCL
#   include "pcl/point_cloud.h"
#   include "pcl/sample_consensus/sac_model_line.h"
#   include "pcl/PointIndices.h"
#endif

namespace am
{
    class HoughLine
    {
        public:
            // CONCEPT
            enum { Dim = 2 };
            typedef float Scalar;
            typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;

            // SET, GET
            inline Eigen::Map<const VectorType> pos() const { return Eigen::Map< const VectorType >( _r_and_angle.data(), Dim ); }
            inline Scalar      & R    ()       { return _r_and_angle(0); }
            inline Scalar const& R    () const { return _r_and_angle(0); }
            inline Scalar      & Angle()       { return _r_and_angle(1); }
            inline Scalar const& Angle() const { return _r_and_angle(1); }

            // CONSTRUCT
            HoughLine()
                : _r_and_angle( Eigen::Matrix<Scalar,Dim,1>::Zero() ) {}
            HoughLine( Scalar r, Scalar angle )
                : _r_and_angle( (VectorType()<< r, angle).finished() ) {}

            template <typename T> static HoughLine
            from3DSpace( std::vector<T> const& coeffs );

            template <typename T>
            static Eigen::Matrix<T,6,1>
            to3DSpace( Eigen::Matrix<T,Eigen::Dynamic,1> r_and_angle );

            template <typename T>
            static std::vector<T>
            to3DSpaceVector( Eigen::Matrix<T,Eigen::Dynamic,1> r_and_angle );

            template <typename T>
            std::vector<T>
            to3DSpaceVector() { return to3DSpaceVector<T>( _r_and_angle ); }

        protected:
            VectorType _r_and_angle;
    };

#ifndef __GF2_HOUGHLINE_HPP__
#define __GF2_HOUGHLINE_HPP__

    template <typename T> HoughLine
    HoughLine::from3DSpace( std::vector<T> const& coeffs )
    {
        Eigen::Vector3f x0( 0.f, 0.f, coeffs[2] );
        Eigen::Vector3f p0 ( coeffs[0], coeffs[1], coeffs[2] );
        Eigen::Vector3f dir( coeffs[3], coeffs[4], coeffs[5] );
        dir.normalize();

        Scalar r     = (x0-p0).cross( dir ).norm(); // / dir.norm() (== 1.f);
        Scalar angle = atan2( dir(1), dir(0) );
        if ( angle < M_PI_4 )   angle += M_PI;

        return HoughLine( r, angle );
    }

    template <typename T>
    std::vector<T>
    HoughLine::to3DSpaceVector( Eigen::Matrix<T,Eigen::Dynamic,1> r_and_angle )
    {
        Eigen::Matrix<Scalar,6,1> eigen_line = to3DSpace( r_and_angle );

        std::vector<Scalar> line_coeff( 6, 0 );
        line_coeff[0] = eigen_line(0);
        line_coeff[1] = eigen_line(1);
        line_coeff[2] = eigen_line(2);
        line_coeff[3] = eigen_line(0);
        line_coeff[4] = eigen_line(1);
        line_coeff[5] = eigen_line(2);

        return line_coeff;
    }

    template <typename T>
    Eigen::Matrix<T,6,1>
    HoughLine::to3DSpace(Eigen::Matrix<T,Eigen::Dynamic,1> r_and_angle )
    {
        Eigen::Matrix<Scalar,6,1> line_coeff;

        Eigen::Matrix<Scalar,3,1> p0( Eigen::Matrix<Scalar,3,1>::Zero() ),
                                  p1( Eigen::Matrix<Scalar,3,1>::Zero() );

        Scalar sin_ang = sin( r_and_angle(1) );

        // debug
        if ( sin_ang != sin_ang )
        {
            std::cerr << "sin for " << r_and_angle(1) << " NAN " << std::endl;
            sin_ang = static_cast<Scalar>( 0 );
        }

        Scalar cos_ang = cos( r_and_angle(1) );

        // debug
        if ( cos_ang != cos_ang )
        {
            std::cerr << "cos for " << r_and_angle(1) << " NAN " << std::endl;
            cos_ang = static_cast<Scalar>( 0 ); // prevent nan
        }

        p0(0) =        (r_and_angle(0) - p0(1)) / sin_ang; // y-crossing
        p1(1) = -1.f * (r_and_angle(0) - p1(0)) / cos_ang; // x crossing

        // debug
        if ( p0(0) != p0(0) )
        {
            std::cerr << "nan, because sin_ang == " << sin_ang << ", r = " << r_and_angle(0) << ", ang = " << r_and_angle(1) << std::endl;
        }
        if ( p1(1) != p1(1) )
        {
            std::cerr << "nan, because cos_ang == " << cos_ang << ", r = " << r_and_angle(0) << ", ang = " << r_and_angle(1) << std::endl;
        }

        //Eigen::Vector3f dir(p1 - p0); dir.normalize();
        line_coeff.template segment<3>(3) = (p1 - p0).normalized();

        // get the closer intersection
        if ( p0.norm() > p1.norm() )
            std::swap(p0,p1);

        line_coeff.template segment<3>(0) = p0;

        return line_coeff;
    }

#endif //__GF2_HOUGHLINE_HPP__

} // ns am

#endif // __GF2_HOUGHLINE_H__

