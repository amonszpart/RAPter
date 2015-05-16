#ifndef RAPTER_ANGLE_HPP
#define RAPTER_ANGLE_HPP

namespace rapter {

    //template<typename Scalar, int Dim> inline Scalar
    //angleInRad( Eigen::Matrix<Scalar,Dim,1> const& v1, Eigen::Matrix<Scalar,Dim,1> const& v2 )
    template<typename Derived, typename DerivedB> inline typename Derived::Scalar
    angleInRad( Derived const& v1, DerivedB const& v2 )
    {
        typedef typename DerivedB::Scalar Scalar;
        // better, than acos
        Scalar angle = atan2( v1.cross(v2).norm(), v1.dot(v2) );

        // fix nans
        if ( angle != angle )   angle = Scalar(0);

        return angle;
    }
    #if 0
    template<typename Scalar, int Dim> inline Scalar
    angleInRadSigned( Eigen::Matrix<Scalar,Dim,1> const& v1, Eigen::Matrix<Scalar,Dim,1> const& v2 )
    {
        // better, than acos
        Scalar angle = atan2( v1.cross(v2).norm(), v1.dot(v2) );

        // fix nans
        if ( angle != angle )   angle = Scalar(0);

        Scalar sign = Eigen::Matrix<Scalar,3,1>::UnitZ()

        return angle;
    }
    #endif

} //...ns rapter

#endif // RAPTER_ANGLE_HPP
