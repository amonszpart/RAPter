#ifndef GF2_UTIL_HPP
#define GF2_UTIL_HPP

#include <vector>
#include "Eigen/Dense"
#include "opencv2/core/core.hpp"

namespace GF2 {
namespace util {

template <typename Scalar>
inline ::cv::Point3f
hsv2rgb( ::cv::Point3_<Scalar> const& in )
{
    double      hh, p, q, t, ff;
    long        i;
    ::cv::Point3f out;

    if ( in.y <= 0.f ) // < is bogus, just shuts up warnings
    {
        out.x = in.z;
        out.y = in.z;
        out.z = in.z;
        return out;
    }

    hh = in.x;
    if ( hh >= 360.f ) hh = 0.f;
    hh /= 60.f;
    i = (long)hh;
    ff = hh - i;
    p = in.z * (1.f - in.y);
    q = in.z * (1.f - (in.y * ff));
    t = in.z * (1.f - (in.y * (1.f - ff)));

    switch(i) {
        case 0:
            out.x = in.z;
            out.y = t;
            out.z = p;
            break;
        case 1:
            out.x = q;
            out.y = in.z;
            out.z = p;
            break;
        case 2:
            out.x = p;
            out.y = in.z;
            out.z = t;
            break;

        case 3:
            out.x = p;
            out.y = q;
            out.z = in.z;
            break;
        case 4:
            out.x = t;
            out.y = p;
            out.z = in.z;
            break;

        case 5:
        default:
            out.x = in.z;
            out.y = p;
            out.z = q;
            break;
    }

    return out;
}

// generates n different colours
// assumes hue [0, 360), saturation [0, 100), lightness [0, 100)
template <typename Scalar>
inline std::vector< ::cv::Point3f >
nColoursCv(int n, Scalar scale, bool random_shuffle )
{
    std::vector< ::cv::Point3f > out;
    //srand(time(NULL));

    float step = 360. / n;
    for ( int i = 0; i < n; ++i )
    {
        ::cv::Point3f c;
        c.x = i * step; // hue
        c.y = (90.f + rand()/(float)RAND_MAX * 10.f) / 100.f; // saturation
        c.z = (50.f + rand()/(float)RAND_MAX * 50.f) / 100.f; // value

        // convert and store
        out.push_back( hsv2rgb(c) * scale );
    }

    if ( random_shuffle )
        std::random_shuffle( out.begin(), out.end() );

    return out;
}

template <typename Scalar>
inline std::vector< ::Eigen::Vector3f >
nColoursEigen( int n, Scalar scale, bool random_shuffle )
{
    typedef std::vector< ::cv::Point3f > CvPointsT;
    CvPointsT colours = nColoursCv( n, scale, random_shuffle );
    std::vector< ::Eigen::Vector3f > colours_eigen;
    //for ( CvPointsT::value_type const& colour : colours )
    for ( CvPointsT::const_iterator it = colours.begin(); it != colours.end(); ++it )
    {
        //colours_eigen.push_back( (Eigen::Vector3f() << colour.x, colour.y, colour.z).finished() );
        colours_eigen.push_back( (Eigen::Vector3f() << it->x, it->y, it->z).finished() );
    }

    return colours_eigen;
}

inline std::string timestamp2Str()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time ( &rawtime );
    timeinfo = localtime (&rawtime);

    strftime ( buffer, 80, "_%Y%m%d_%H%M", timeinfo );

    return std::string( buffer );
} //...timestamp2Str()

} //...namespace util
} //...namespace GF2

#endif // GF2_UTIL_HPP
