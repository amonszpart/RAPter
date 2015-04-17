#ifndef GO_TRIANGLE_HPP
#define GO_TRIANGLE_HPP

#include "Eigen/Dense"
#include <iostream>

namespace globopt
{

    template <typename _Scalar>
    class Triangle
    {
        public:

            typedef          _Scalar                   Scalar;
            typedef typename Eigen::Matrix<Scalar,3,1> VectorT;
            typedef typename Eigen::Matrix<Scalar,3,3, Eigen::ColMajor> CornersT;

        public:

            inline Triangle() : corners( CornersT::Zero() ) {}
            inline Triangle( VectorT a, VectorT b, VectorT c ) : corners( (CornersT() << a,b,c).finished() ) {}
            // todo: figure out non-copy return type
            inline const VectorT getCorner( int const id ) const { return corners.col( id ); }
            inline Scalar getDistance( Eigen::Matrix<Scalar,3,1> const& point ) const { return getDistance(*this,point); }
            inline Scalar getSquaredDistance( Eigen::Matrix<Scalar,3,1> const& point ) const { return getSquaredDistance(*this,point); }
            template <typename _VectorT>
            static inline Scalar getSquaredDistance( Triangle const& tri, _VectorT const& point );
            template <typename _VectorT>
            static inline Scalar getDistance( Triangle const& tri, _VectorT const& point );
            Triangle<_Scalar>::VectorT dir() const;
            VectorT getMean() const { return corners.rowwise().mean(); }

        protected:
            //! \brief Corner points in columns
            CornersT corners;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    }; //...cls Triangle

} //...ns globopt

namespace globopt
{
    template <typename _Scalar>
    inline typename Triangle<_Scalar>::VectorT Triangle<_Scalar>::dir() const
    {
        VectorT edge0 = this->getCorner(1) - this->getCorner(0);
//        std::cout << "edge0: " << edge0.template transpose()  << " = "
//                              << this->getCorner(1).template transpose() << " - "
//                              << this->getCorner(0).template transpose()
//                              << std::endl;
        VectorT edge1 = this->getCorner(2) - this->getCorner(0);
//        std::cout << "edge1: " << edge1.template transpose()  << " = " << this->getCorner(2).template transpose() << " - " << this->getCorner(0).template transpose()
//                              << std::endl;
        //return edge0.normalized().cross( edge1.normalized() );
        return edge0.cross( edge1 ).normalized();
    } //...dir()

    // from: http://kyousai.googlecode.com/svn/trunk/LibMathematics/Distance/Wm5DistPoint3Triangle3.cpp
    template <typename _Scalar>
    template <typename _VectorT>
    inline _Scalar Triangle<_Scalar>::getSquaredDistance( Triangle<_Scalar> const& tri, _VectorT const& point )
    {

        _VectorT diff  = tri.getCorner(0) - point;
        _VectorT edge0 = tri.getCorner(1) - tri.getCorner(0);
        _VectorT edge1 = tri.getCorner(2) - tri.getCorner(0);
        Scalar a00 = edge0.squaredNorm();
        Scalar a01 = edge0.dot(edge1);
        Scalar a11 = edge1.squaredNorm();
        Scalar b0 = diff.dot(edge0);
        Scalar b1 = diff.dot(edge1);
        Scalar c = diff.squaredNorm();
        Scalar det = std::abs(a00*a11 - a01*a01);
        Scalar s = a01*b1 - a11*b0;
        Scalar t = a01*b0 - a00*b1;
        Scalar sqrDistance;

        if (s + t <= det)
        {
            if (s < (Scalar)0)
            {
                if (t < (Scalar)0)  // region 4
                {
                    if (b0 < (Scalar)0)
                    {
                        t = (Scalar)0;
                        if (-b0 >= a00)
                        {
                            s = (Scalar)1;
                            sqrDistance = a00 + ((Scalar)2)*b0 + c;
                        }
                        else
                        {
                            s = -b0/a00;
                            sqrDistance = b0*s + c;
                        }
                    }
                    else
                    {
                        s = (Scalar)0;
                        if (b1 >= (Scalar)0)
                        {
                            t = (Scalar)0;
                            sqrDistance = c;
                        }
                        else if (-b1 >= a11)
                        {
                            t = (Scalar)1;
                            sqrDistance = a11 + ((Scalar)2)*b1 + c;
                        }
                        else
                        {
                            t = -b1/a11;
                            sqrDistance = b1*t + c;
                        }
                    }
                }
                else  // region 3
                {
                    s = (Scalar)0;
                    if (b1 >= (Scalar)0)
                    {
                        t = (Scalar)0;
                        sqrDistance = c;
                    }
                    else if (-b1 >= a11)
                    {
                        t = (Scalar)1;
                        sqrDistance = a11 + ((Scalar)2)*b1 + c;
                    }
                    else
                    {
                        t = -b1/a11;
                        sqrDistance = b1*t + c;
                    }
                }
            }
            else if (t < (Scalar)0)  // region 5
            {
                t = (Scalar)0;
                if (b0 >= (Scalar)0)
                {
                    s = (Scalar)0;
                    sqrDistance = c;
                }
                else if (-b0 >= a00)
                {
                    s = (Scalar)1;
                    sqrDistance = a00 + ((Scalar)2)*b0 + c;
                }
                else
                {
                    s = -b0/a00;
                    sqrDistance = b0*s + c;
                }
            }
            else  // region 0
            {
                // minimum at interior point
                Scalar invDet = ((Scalar)1)/det;
                s *= invDet;
                t *= invDet;
                sqrDistance = s*(a00*s + a01*t + ((Scalar)2)*b0) +
                    t*(a01*s + a11*t + ((Scalar)2)*b1) + c;
            }
        }
        else
        {
            Scalar tmp0, tmp1, numer, denom;

            if (s < (Scalar)0)  // region 2
            {
                tmp0 = a01 + b0;
                tmp1 = a11 + b1;
                if (tmp1 > tmp0)
                {
                    numer = tmp1 - tmp0;
                    denom = a00 - ((Scalar)2)*a01 + a11;
                    if (numer >= denom)
                    {
                        s = (Scalar)1;
                        t = (Scalar)0;
                        sqrDistance = a00 + ((Scalar)2)*b0 + c;
                    }
                    else
                    {
                        s = numer/denom;
                        t = (Scalar)1 - s;
                        sqrDistance = s*(a00*s + a01*t + ((Scalar)2)*b0) +
                            t*(a01*s + a11*t + ((Scalar)2)*b1) + c;
                    }
                }
                else
                {
                    s = (Scalar)0;
                    if (tmp1 <= (Scalar)0)
                    {
                        t = (Scalar)1;
                        sqrDistance = a11 + ((Scalar)2)*b1 + c;
                    }
                    else if (b1 >= (Scalar)0)
                    {
                        t = (Scalar)0;
                        sqrDistance = c;
                    }
                    else
                    {
                        t = -b1/a11;
                        sqrDistance = b1*t + c;
                    }
                }
            }
            else if (t < (Scalar)0)  // region 6
            {
                tmp0 = a01 + b1;
                tmp1 = a00 + b0;
                if (tmp1 > tmp0)
                {
                    numer = tmp1 - tmp0;
                    denom = a00 - ((Scalar)2)*a01 + a11;
                    if (numer >= denom)
                    {
                        t = (Scalar)1;
                        s = (Scalar)0;
                        sqrDistance = a11 + ((Scalar)2)*b1 + c;
                    }
                    else
                    {
                        t = numer/denom;
                        s = (Scalar)1 - t;
                        sqrDistance = s*(a00*s + a01*t + ((Scalar)2)*b0) +
                            t*(a01*s + a11*t + ((Scalar)2)*b1) + c;
                    }
                }
                else
                {
                    t = (Scalar)0;
                    if (tmp1 <= (Scalar)0)
                    {
                        s = (Scalar)1;
                        sqrDistance = a00 + ((Scalar)2)*b0 + c;
                    }
                    else if (b0 >= (Scalar)0)
                    {
                        s = (Scalar)0;
                        sqrDistance = c;
                    }
                    else
                    {
                        s = -b0/a00;
                        sqrDistance = b0*s + c;
                    }
                }
            }
            else  // region 1
            {
                numer = a11 + b1 - a01 - b0;
                if (numer <= (Scalar)0)
                {
                    s = (Scalar)0;
                    t = (Scalar)1;
                    sqrDistance = a11 + ((Scalar)2)*b1 + c;
                }
                else
                {
                    denom = a00 - ((Scalar)2)*a01 + a11;
                    if (numer >= denom)
                    {
                        s = (Scalar)1;
                        t = (Scalar)0;
                        sqrDistance = a00 + ((Scalar)2)*b0 + c;
                    }
                    else
                    {
                        s = numer/denom;
                        t = (Scalar)1 - s;
                        sqrDistance = s*(a00*s + a01*t + ((Scalar)2)*b0) +
                            t*(a01*s + a11*t + ((Scalar)2)*b1) + c;
                    }
                }
            }
        }

        // Account for numerical round-off error.
        if (sqrDistance < (Scalar)0)
        {
            sqrDistance = (Scalar)0;
        }

    //        mClosestPoint0 = point;
    //        mClosestPoint1 = tri.getCorner(0) + s*edge0 + t*edge1;
    //        mTriangleBary[1] = s;
    //        mTriangleBary[2] = t;
    //        mTriangleBary[0] = (Scalar)1 - s - t;
        return sqrDistance;
    } //...getSquaredDistance()

    template <typename _Scalar>
    template <typename _VectorT>
    inline _Scalar Triangle<_Scalar>::getDistance( Triangle<_Scalar> const& tri, _VectorT const& point ) { return std::sqrt( Triangle<_Scalar>::getSquaredDistance(tri,point) ); }
}


#if 0

int testTriangle()
{
    typedef typename Eigen::Matrix<Scalar,3,1> Vector;
    typedef Triangle<Scalar> Triangle;

    std::vector< std::pair<Vector, Scalar> > testPoints = { {Vector(0.5,0.5,0), 0.},
                                                            {Vector(0.,0.,0.), 0.},
                                                            {Vector(0.,1.,0.), 0.},
                                                            {Vector(1.,0.,0.), 0.},
                                                            {Vector(1.,0.,1.), 1.},
                                                            {Vector(-1.,-1,-1.), std::sqrt(3.)},
                                                            {Vector(0.45,-1,-1.), std::sqrt(2.)}
                                                          };
    Triangle tri( Vector(0,0,0), Vector(0,1,0), Vector(1,0,0) );
    for ( auto it = testPoints.begin(); it != testPoints.end(); ++it )
    {
        Scalar dist = getDistance<Scalar>( tri, it->first );
        Scalar diff = std::abs( dist - it->second );
        if ( diff < Scalar(0.01) )
            std::cout << dist << " == " << it->second << std::endl;
        else
            std::cerr << dist << " == " << it->second << std::endl;
    }
}

return EXIT_FAILURE;
#endif

#endif // GO_TRIANGLE_HPP
