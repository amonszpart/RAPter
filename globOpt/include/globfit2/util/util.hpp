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

static void RGB2HSV(float r, float g, float b,
                    float &h, float &s, float &v)
{
    float K = 0.f;

    if (g < b)
    {
        std::swap(g, b);
        K = -1.f;
    }

    if (r < g)
    {
        std::swap(r, g);
        K = -2.f / 6.f - K;
    }

    float chroma = r - std::min(g, b);
    h = fabs(K + (g - b) / (6.f * chroma + 1e-20f));
    s = chroma / (r + 1e-20f);
    v = r;
}

/*! \brief
 *  \tparam MatrixDervied Concept: Eigen::Matrix<Scalar,3,1>
 */
template <class MatrixDerived>
inline MatrixDerived
hsv2rgbEigen( MatrixDerived const& in )
{
    typedef typename MatrixDerived::Scalar Scalar;

    double      hh, p, q, t, ff;
    long        i;
    MatrixDerived out;

    if ( in(1) <= Scalar(0.) ) // < is bogus, just shuts up warnings
    {
        out(0) = in(2);
        out(1) = in(2);
        out(2) = in(2);
        return out;
    }

    hh = in(0);
    if ( hh >= Scalar(360.) ) hh = Scalar(0.);
    hh /= Scalar(60.);
    i = (long)hh;
    ff = hh - i;
    p = in(2) * (1.f - in(1));
    q = in(2) * (1.f - (in(1) * ff));
    t = in(2) * (1.f - (in(1) * (1.f - ff)));

    switch(i) {
        case 0:
            out(0) = in(2);
            out(1) = t;
            out(2) = p;
            break;
        case 1:
            out(0) = q;
            out(1) = in(2);
            out(2) = p;
            break;
        case 2:
            out(0) = p;
            out(1) = in(2);
            out(2) = t;
            break;

        case 3:
            out(0) = p;
            out(1) = q;
            out(2) = in(2);
            break;
        case 4:
            out(0) = t;
            out(1) = p;
            out(2) = in(2);
            break;

        case 5:
        default:
            out(0) = in(2);
            out(1) = p;
            out(2) = q;
            break;
    }

    return out;
}

template <typename Scalar> inline
static void rgb2hsv( Scalar r, Scalar g, Scalar b,
                     Scalar &h, Scalar &s, Scalar &v)
{
    Scalar K(0.);
    if (g < b)
    {
        std::swap(g, b);
        K = -1.;
    }
    Scalar min_gb(b);
    if (r < g)
    {
        std::swap(r, g);
        K = -2. / 6. - K;
        min_gb = std::min(g, b);
    }
    Scalar chroma = r - min_gb;
    h = fabs(K + (g - b) / (6. * chroma + 1.e-20));
    s = chroma / (r + 1.e-20);
    v = r;
}

// generates n different colours
// assumes hue [0, 360), saturation [0, 100), lightness [0, 100)
template <typename Scalar>
inline std::vector< ::cv::Point3f >
nColoursCv(int n, Scalar scale, bool random_shuffle, float min_value = 50.f, float min_saturation = 90.f )
{
    std::vector< ::cv::Point3f > out;
    //srand(time(NULL));

    float step = 360. / n;
    const float saturation_rest = 100.f - min_saturation;
    const float value_rest = 100.f - min_value;
    for ( int i = 0; i < n; ++i )
    {
        ::cv::Point3f c;
        c.x = i * step; // hue
        c.y = (min_saturation + rand()/(float)RAND_MAX * saturation_rest) / 100.f; // saturation
        c.z = (min_value      + rand()/(float)RAND_MAX * value_rest     ) / 100.f; // value

        // convert and store
        out.push_back( hsv2rgb(c) * scale );
    }

    if ( random_shuffle )
        std::random_shuffle( out.begin(), out.end() );

    return out;
}

template <typename Scalar>
inline std::vector< ::Eigen::Vector3f >
nColoursEigen( int n, Scalar scale, bool random_shuffle, float min_value = 50.f, float min_saturation = 90.f )
{
    typedef std::vector< ::cv::Point3f > CvPointsT;
    CvPointsT colours = nColoursCv( n, scale, random_shuffle, min_value, min_saturation );
    std::vector< ::Eigen::Vector3f > colours_eigen;
    //for ( CvPointsT::value_type const& colour : colours )
    for ( CvPointsT::const_iterator it = colours.begin(); it != colours.end(); ++it )
    {
        //colours_eigen.push_back( (Eigen::Vector3f() << colour.x, colour.y, colour.z).finished() );
        colours_eigen.push_back( (Eigen::Vector3f() << it->x, it->y, it->z).finished() );
    }

    return colours_eigen;
}

inline std::vector< ::Eigen::Vector3f >
paletteLightColoursEigen( int min_count = 0 )
{
    // don't have template parameters, since we return Vector3f....
    std::vector< ::Eigen::Vector3f > colours_eigen;
    colours_eigen.resize(7);
    colours_eigen [2] << 184.f, 210.f, 236.f;
    colours_eigen [1] << 217.f, 228.f, 170.f;
    colours_eigen [0] << 242.f, 175.f, 173.f;
    colours_eigen [3] << 243.f, 209.f, 176.f;
    colours_eigen [4] << 213.f, 178.f, 212.f;
    colours_eigen [5] << 221.f, 185.f, 169.f;
    colours_eigen [6] << 235.f, 192.f, 218.f;

    std::vector< ::Eigen::Vector3f > out;
    while ( out.size() < min_count )
        out.insert( out.end(), colours_eigen.begin(), colours_eigen.end() );

    return out;
}

inline Eigen::Vector3f
paletteLightNeutralColour(){
    return Eigen::Vector3f(204.f, 204.f, 204.f);
}

inline std::vector< ::Eigen::Vector3f >
paletteMediumColoursEigen( int min_count = 0 )
{
    // don't have template parameters, since we return Vector3f....
    std::vector< ::Eigen::Vector3f > colours_eigen;
    colours_eigen.resize(7);
    colours_eigen [2] << 090.f, 155.f, 212.f;
    colours_eigen [1] << 122.f, 195.f, 106.f;
    colours_eigen [0] << 241.f, 090.f, 096.f;
    colours_eigen [3] << 250.f, 167.f, 091.f;
    colours_eigen [4] << 158.f, 103.f, 171.f;
    colours_eigen [5] << 206.f, 112.f, 088.f;
    colours_eigen [6] << 215.f, 127.f, 180.f;

    std::vector< ::Eigen::Vector3f > out;
    while ( out.size() < min_count )
        out.insert( out.end(), colours_eigen.begin(), colours_eigen.end() );

    return out;
}

inline Eigen::Vector3f
paletteMediumNeutralColour(){
    return Eigen::Vector3f(115.f, 115.f, 115.f);
}

inline std::vector< ::Eigen::Vector3f >
paletteDarkColoursEigen( int min_count = 0 )
{
    // don't have template parameters, since we return Vector3f....
    std::vector< ::Eigen::Vector3f > colours_eigen;
    colours_eigen.resize(7);
    colours_eigen [2] << 024.f, 090.f, 169.f;
    colours_eigen [1] << 000.f, 140.f, 072.f;
    colours_eigen [0] << 238.f, 046.f, 047.f;
    colours_eigen [3] << 244.f, 125.f, 035.f;
    colours_eigen [4] << 102.f, 044.f, 145.f;
    colours_eigen [5] << 162.f, 029.f, 033.f;
    colours_eigen [6] << 180.f, 056.f, 148.f;

    std::vector< ::Eigen::Vector3f > out;
    while ( out.size() < min_count )
        out.insert( out.end(), colours_eigen.begin(), colours_eigen.end() );

    return out;
}

inline Eigen::Vector3f
paletteDarkNeutralColour(){
    return Eigen::Vector3f(001.f, 002.f, 002.f);
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

inline int parseIteration( std::string const& input_prims_path )
{
    int iteration = 0;

    size_t it_loc = input_prims_path.find("_it");
    if ( it_loc == std::string::npos )
        iteration = -1;
    else
    {
        iteration = atoi( input_prims_path.substr( it_loc+3,
                                                       input_prims_path.find(".",it_loc)).c_str() );
        //iteration = atoi( input_prims_path.substr( it_loc+3,1 ).c_str() );
    }

    return iteration;
}

#if 0
template <class ExtremaT, class PrimitiveT>
void printme (ExtremaT& extre, PrimitiveT& prim, std::ofstream& mstream)
{
    using std::endl;

    mstream << extre.rbegin()->transpose() << endl;
    mstream << extre.begin()->transpose() << endl << endl << endl;

    auto mit_local = extre.begin();
    auto mit_localn = extre.begin()+1;

    for ( ; mit_localn != extre.end(); ++mit_local, ++mit_localn )
    {
        mstream << mit_local->transpose() << endl;
        mstream << mit_localn->transpose() << endl << endl << endl;
    }
    typedef typename ExtremaT::value_type PointT;
    typedef typename PointT::Scalar _Scalar;

    PointT center (PointT::Zero());
    std::for_each(extre.begin(), extre.end(), [&center] (const PointT& p){ center+=p; });
    center /= _Scalar(extre.size());

    // display normal
    mstream << center.transpose() << endl;
    mstream << (center+prim.normal()).transpose() << endl << endl << endl;

    // display normal from extrema
    mstream << center.transpose() << endl;
    mstream << (center+(extre[1]-extre[0]).normalized().cross((extre[2]-extre[1])).normalized()).transpose() << endl << endl << endl;
}

std::ofstream f("onplanecloud.plot");
for ( int i = 0; i != on_plane_cloud.size(); ++i )
{
    f << on_plane_cloud[i].template pos().transpose() << "\n";
}
f.close();
std::ofstream f2("onplanecloud_normal.plot");
f2 << this->pos().transpose() << "\n"
   << (this->pos() + this->dir()).transpose() << "\n";
f2 << "\n\n"
   << this->pos().transpose() << "\n"
   << (this->pos() + (minMax[2] - minMax[1]).normalized().cross( (minMax[1] - minMax[0]).normalized() ).normalized()).transpose() << "\n";

f2.close();
system( "gnuplot -e \"splot 'onplanecloud.plot' w points, 'onplanecloud_normal.plot' w lines\" -p" );
#endif

} //...namespace util
} //...namespace GF2

#endif // GF2_UTIL_HPP
