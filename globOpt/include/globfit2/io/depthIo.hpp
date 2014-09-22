#ifndef GF2_DEPTHIO_HPP
#define GF2_DEPTHIO_HPP

#include <string>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include "Eigen/Dense"

namespace GF2 {
namespace io {

/*! \brief Load a depth image stored as a .dat file.
 *         Assumes, that first two int entries are height and width (in this order).
 * \tparam     FilePixelT The scalar format of the data file Concept: float.
 * \tparam     MatPixelT  The output format of cv::Mat. Concept: ushort.
 * \param[out] mat        Output opencv matrix.
 * \param[in]  path       File path to read from.
 * \param[in]  alpha      Multiplier to apply to each read pixel.
*/
template <typename MatPixelT, typename FilePixelT> inline int
loadDepthAsBinary( ::cv::Mat &mat, std::string path, float alpha = 1.f )
{
    FILE *fp = fopen( path.c_str(), "rb" );
    if (!fp)
    {
        printf( "loadDepthAsBinary failed, when opening %s\n", path.c_str() );
        return false;
    }

    int w,h, read_size = 0;
    read_size += fread( &h,sizeof(int),1,fp );
    read_size += fread( &w,sizeof(int),1,fp );
    if ( read_size != 2 )
    {
        std::cerr << "[" << __func__ << "]: " << "w or h not read properly: " << w << ", " << h << std::endl;
    }

    FilePixelT *p_depth = new FilePixelT[ h * w ];

    read_size = fread( p_depth, sizeof(FilePixelT), w*h, fp );
    if ( read_size != w*h )
        std::cerr << "[" << __func__ << "]: " << "read size != " << w*h << std::endl;

    mat.create( h, w, CV_16UC1 );
    int p_depth_index = 0;
    for ( unsigned y = 0; y < mat.rows; ++y  )
    {
        for ( unsigned x = 0; x < mat.cols; ++x, ++p_depth_index )
        {
            mat.at<MatPixelT>(y,x) = p_depth[ p_depth_index ] * alpha;
        }
    }

    fclose( fp );
    if ( p_depth )
        delete [] p_depth;

    return EXIT_SUCCESS;
} //...loadDepthAsBinary()

/*! \brief Loads a raw, dat or png depth map, and outputs [mm] values.
 *
 *  \param[in] depth_path File path to read depth map from.
 *  \return OpenCV 2D matrix with ushort depth values in [mm].
 */
cv::Mat
loadDepth( std::string depth_path )
{
    cv::Mat dep;
    if ( depth_path.find("raw") != std::string::npos )
    {
        loadDepthAsBinary<ushort,ushort>( dep, depth_path, 1. );
    }
    else if ( depth_path.find("dat") != std::string::npos )
    {
        // NYU, read float, and convert it to ushort
        loadDepthAsBinary<ushort,float>( dep, depth_path, 1000. );
    }
    else
    {
        dep = cv::imread( depth_path, cv::IMREAD_UNCHANGED );
    }

    return dep;
} //...loadDepth()

// Aron's kinect calibration
#define FX_D 588.486898
#define FY_D 588.729416
#define CX_D 314.509516
#define CY_D 242.997237
#define K1_D -0.094563
#define K2_D 0.387595
#define P1_D 0.000838
#define P2_D 0.001505
#define K3_D -0.494912
#define ALPHA_D -0.002096

template <typename _Scalar>
struct Intrinsics
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef _Scalar Scalar;
    enum CALIB_ID { TIANJIAS, ARONS, /*ARONS_WDISTORTION,*/ RGBDEMOS };

    inline Intrinsics( CALIB_ID calib_id = ARONS )
    {
         switch ( calib_id )
         {
             case ARONS:
                 _intrinsics << FX_D, 0., CX_D,
                                0., FY_D, CY_D,
                                0., 0., 1.;
                 break;
             case TIANJIAS:
                 _intrinsics << 5.1930334103339817e+02, 0., 3.2850951551345941e+02,
                                0., 5.1816401430246583e+02, 2.5282555217253503e+02,
                                0., 0., 1.;
                 break;
             case RGBDEMOS:
                 _intrinsics << 5.4013644168716110e+02, 0., 320.,
                                0., 5.4013644168716110e+02, 240.,
                                0., 0., 1.;
                 break;
         }
    }

    inline Intrinsics( Scalar fx, Scalar fy, Scalar cx, Scalar cy )
        : _intrinsics( (Eigen::Matrix<Scalar,3,3>() <<
                        fx, 0., cx,
                        0., fy, cy,
                        0., 0., 1.).finished() ) {}
    operator Eigen::Matrix<_Scalar,3,3>() { return _intrinsics; }
    protected:
        Eigen::Matrix<Scalar,3,3> _intrinsics;
};

template <typename _Scalar> inline Eigen::Matrix<_Scalar,3,1>
point2To3D( Eigen::Matrix<_Scalar,2,1> const& pnt2, Eigen::Matrix<_Scalar,3,3> const& intrinsics )
{
    return (Eigen::Matrix<_Scalar,3,1>() << (pnt2(0) - intrinsics(0,2)) / intrinsics(0,0), // (x - cx) / fx
                                            (pnt2(1) - intrinsics(1,2)) / intrinsics(1,1), // (y - cy) / fy
                                            1.f ).finished();
}

/*! \brief Decides, if a depth map pixel is a valid 3D point.
 *  \return true, if depth > 20cm. (> 0.2)
 */
template <typename _Scalar>
inline bool isValidDepth( _Scalar d )
{
    if ( d > _Scalar(10.) ) std::cerr << "[" << __func__ << "]: " << "Depth values are assumed to be between 0..10, and not " << d << std::endl;

    return d > _Scalar(.2);
}

/*! \brief Merges a depth image and an optional rgb image to a point cloud.
 */
template <typename depT, typename _Scalar, typename _PclCloudT>
inline int rgbd2PointCloud( _PclCloudT                                     & pcl_cloud
                          , cv::Mat                                   const& dep
                          , cv::Mat                                   const& img
                          , _Scalar                                   const  alpha      = 1. / 1000.
                          , Eigen::Matrix<_Scalar,3,3>                const& intrinsics = Intrinsics<_Scalar>()
                          , Eigen::Transform<_Scalar,3,Eigen::Affine> const* pose       = NULL
                          )
{
    typedef Eigen::Matrix<_Scalar,3,1>            Position;
    typedef typename _PclCloudT::PointType        PclPointT;
    //typedef typename _PointContainerT::value_type PointPrimitiveT;

    // check inputs to be the same size
    if ( !img.empty() )
    {
        if ( dep.size() != img.size() )
        {
            std::cerr << "matsTo3D(): dep and rgb need the same size! "
                      << dep.rows << "x" << dep.cols << ", "
                      << img.rows << "x" << img.cols << std::endl;
            return EXIT_FAILURE;
        }
    }

    // check image dimensions
    if ( (!img.empty()) && (img.channels() != 3) )
    {
        std::cerr << "matsTo3D(): rgb should be 3 channels! " << img.channels() << std::endl;
        return EXIT_FAILURE;
    }

    // allocate ouptut
    //points = PCloud::Ptr( new PCloud );
//    pcl_cloud.clear();
//    pcl_cloud.reserve( img.cols * img.rows );

    // point to read to, and insert
    Position pnt3D;
    uint32_t rgb = 128 << 16 | 128 << 8 | 128;

    Eigen::Matrix<_Scalar,3,3> rotation;
    Eigen::Matrix<_Scalar,3,1> translation;
    if ( pose )
    {
        rotation    = pose->rotation();
        translation = pose->translation();
    }

    // copy inputs
    for ( int y = 0; y < dep.rows; ++y )
    {
        for ( int x = 0; x < dep.cols; ++x )
        {
            //  read depth
            float depth = (_Scalar)dep.at<depT>( y,x ) * alpha;
            // check validity
            bool isValid = isValidDepth( depth );
            if ( !isValid )
                continue;

            // calculate point
            pnt3D = point2To3D( (Eigen::Matrix<_Scalar,2,1>() << x,y).finished(), intrinsics ) * depth;

            // apply pose
            if ( pose )
                pnt3D = rotation * pnt3D + translation;

#if 0
            pcl_cloud.push_back( PointPrimitiveT( /* pos: */ pnt3D, /* dir: */ (Position()<<1.,1.,1.).finished()) );
#else
            // convert
            PclPointT point;
            point.x = pnt3D( 0 );
            point.y = pnt3D( 1 );
            point.z = pnt3D( 2 );

            // add colour
            if ( !img.empty() )
            {
                rgb = (static_cast<uint32_t>(img.at<uchar>(y, x * img.channels() + 2)) << 16 |
                       static_cast<uint32_t>(img.at<uchar>(y, x * img.channels() + 1)) << 8  |
                       static_cast<uint32_t>(img.at<uchar>(y, x * img.channels()    ))        );
            }
            point.rgb = *reinterpret_cast<float*>( &rgb );

            // add point
            pcl_cloud.push_back( point );
#endif
        }
    }

    //points->width = (int)points->points.size ();
    //points->height = 1;

    return EXIT_SUCCESS;

} // matsTo3D

} //...ns io
} //...ns GF2

#endif // GF2_DEPTHIO_HPP
