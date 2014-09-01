#ifndef KMEANS_H
#define KMEANS_H

#include <vector>
#include "my_types.h"
#include "KMPoint.h"
#include "primitive.h"

namespace am
{
    namespace kmeans
    {

        // concept
        class My3DPoint
        {
            public:
                enum { Dim = 3 };
                typedef float Scalar;
                typedef Eigen::Matrix<Scalar,Dim,1> VectorType;

                My3DPoint()
                    : _pnt(0.f,0.f,0.f)
                {}

                My3DPoint( MyPoint const& pnt )
                    : _pnt( pnt )
                {}

                inline Scalar & X() { return _pnt.x; }
                inline Scalar & Y() { return _pnt.y; }
                inline Scalar & Z() { return _pnt.z; }
                inline Scalar const& X() const { return _pnt.x; }
                inline Scalar const& Y() const { return _pnt.y; }
                inline Scalar const& Z() const { return _pnt.z; }

                inline Eigen::Map<const VectorType>  pos() const { return Eigen::Map< const VectorType >( _pnt.data, Dim ); }

            protected:
                MyPoint                _pnt;
        };

        extern double randf(double m);

        /**
         * TPoint::Scalar
         * TPoint::Dim
         * TPoint::VectorType      Eigen::Matrix<Scalar, Dim, 1>
         * TPoint::pos()           Eigen::Map<Scalar,Dim,1> pos()
         * TPoints::const_iterator
         * TPoints::begin()
         * TPoints::end()
         * TPoints::size()
         */
        template <class TPoints, class TPoint>
        class KMeans
        {
            public:
                //enum { PRIMITIVE_POINT = 1, PRIMITIVE_POINT_1D = 2 };
                typedef typename TPoint::VectorType VectorType;
                typedef typename TPoint::Scalar     Scalar;
                typedef typename TPoints::const_iterator         ConstIterator;


                static int
                cluster( typename Primitive<TPoints,TPoint>::Vec & primitives
                         , std::vector<int>                      & labels
                         , TPoints                          const& points
                         , VectorType                       const* dim_coeffs
                         , bool                                    kmeans_pp_primitives = false );

                static int
                toJson( std::string           path
                        , TPoint          const* p_points
                        , int                 points_size
                        , std::vector<int>                  & labels
                        , typename Primitive<TPoints,TPoint>::Vec   const& primitives
                        , std::vector<Eigen::Vector3f> const& colours
                        );

            protected:
                static void
                kpp( typename Primitive<TPoints,TPoint>::Vec & primitives
                     , std::vector<int>                      & labels
                     , TPoints                          const& points
                     , VectorType                       const* dim_coeffs );
                static int
                lloyd( typename Primitive<TPoints,TPoint>::Vec & primitives
                       , std::vector<int>                      & labels
                       , TPoints                          const& points
                       , VectorType                       const* dim_coeffs );
                static int
                nearest_primitive( TPoint                                    const& pt
                                   , int                                            curr_label
                                   , typename Primitive<TPoints,TPoint>::Vec const& primitives
                                   , Scalar                                       * depth_out         = NULL
                                                                                                        , VectorType                              const* dim_coeffs        = NULL
                                                                                                                                                                             , typename Primitive<TPoints,TPoint>::Vec::const_iterator *end_arg = NULL);
        };

    } // ns kmeans

} // ns am

#ifndef INC_KMEANS_HPP
#   define INC_KMEANS_HPP
#   include "kMeans.hpp"
#endif

#endif // KMEANS_H

