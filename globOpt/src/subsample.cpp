#include "Eigen/Dense"
#include "boost/filesystem.hpp"

#include "pcl/console/parse.h"
#include "pcl/io/io.h"
#include "pcl/io/ply_io.h"
#include "pcl/common/common.h"


class CloudColouring
{
    public:
        template <typename Scalar = float> static inline Eigen::Matrix<int,3,1>
        getColourForDistance( Scalar distance, Scalar max_dist, Scalar min_dist );
};

template <typename Scalar> Eigen::Matrix<int,3,1> inline
CloudColouring::getColourForDistance( Scalar distance, Scalar max_dist, Scalar min_dist )
{
    Eigen::Matrix<int,3,1> ret_colour( Eigen::Matrix<int,3,1>::Zero() );

    Scalar spread   = (max_dist - min_dist);            // calculate domain width for normalization
    Scalar val      = (distance - min_dist) / spread;   // unbias, normalize

    if ( (val > Scalar(1)) || (val < Scalar(0)) )
    {
        std::cerr << "distance: " << distance
                  << ", max_dist: " << max_dist
                  << ", min_dist: " << min_dist
                  << ", spread: " << spread
                  << ", val: " << val << std::endl;
    }
    // borrowed from openkinect
    val = 0.22 + val * 0.323996; // map spread of visualization...
          // 0.046007 + val * 0.498993;
    int pval = val * val * val * 36 * 256; //m_gamma[depth[i]];
    int lb   = pval & 0xff;
    switch ( pval>>8 )
    {
        case 0:
            ret_colour(0) = 255;
            ret_colour(1) = 255-lb;
            ret_colour(2) = 255-lb;
            break;
        case 1:
            ret_colour(0) = 255;
            ret_colour(1) = lb;
            ret_colour(2) = 0;
            break;
        case 2:
            ret_colour(0) = 255-lb;
            ret_colour(1) = 255;
            ret_colour(2) = 0;
            break;
        case 3:
            ret_colour(0) = 0;
            ret_colour(1) = 255;
            ret_colour(2) = lb;
            break;
        case 4:
            ret_colour(0) = 0;
            ret_colour(1) = 255-lb;
            ret_colour(2) = 255;
            break;
        case 5:
            ret_colour(0) = 0;
            ret_colour(1) = 0;
            ret_colour(2) = 255-lb;
            break;
        default:
            ret_colour(0) = 0;
            ret_colour(1) = 0;
            ret_colour(2) = 0;
            break;
    }

    return ret_colour;
}

int subsample( int argc, char** argv )
{
    typedef float                     Scalar;
    typedef Eigen::Matrix<Scalar,4,1> Position;

    bool        valid_input             = true;
    std::string cloud_path              = "./cloud.ply";
    float       N                       = 1000.;

    // cloud
    if ( (pcl::console::parse_argument( argc, argv, "--cloud", cloud_path) < 0)
         || !boost::filesystem::exists( cloud_path ) )
    {
        std::cerr << "[" << __func__ << "]: " << "--cloud does not exist: " << cloud_path << std::endl;
        valid_input = false;
    }

    if ( pcl::console::parse_argument( argc, argv, "--N", N) < 0 )
    {
        std::cerr << "[" << __func__ << "]: " << "Need N to work" << std::endl;
        valid_input = false;
    }

    Position sceneSize;
    {
        float scene_size = 1.f;
        if ( pcl::console::parse_argument( argc, argv, "--scene-size", scene_size) < 0 )
        {
            std::cerr << "[" << __func__ << "]: " << "\"--scene-size 1\" assumed to normalize scene to" << std::endl;
        }
        sceneSize << scene_size, scene_size, scene_size, 1;
    }

    Position origin;
    std::vector<Scalar> originCoords;
    bool doOrigin = false;
    if ( pcl::console::parse_x_arguments( argc, argv, "--origin", originCoords ) >= 0 )
    {
        origin << originCoords[0], originCoords[1], originCoords[2], 0;
        doOrigin = true;
    }
    else
        origin << 0,0,0,0;

    if ( !valid_input || (pcl::console::find_switch(argc,argv,"-h")) || (pcl::console::find_switch(argc,argv,"--help")) )
    {
        std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << "--subsample\n"
                  << "\t--cloud " << cloud_path << "\n"
                  << "\t--N " << N
                  << "\t--scene-size " << sceneSize.transpose()
                  << "\t--origin " << origin.transpose() << " \t to colour by distance from origin"
                  << "\n";

        return EXIT_FAILURE;
    }

    srand( 123456 );

    pcl::PointCloud<pcl::PointXYZRGB> cloud, out_cloud;
    out_cloud.reserve( cloud.size() );
    std::cout << "reading " << cloud_path << "..."; fflush(stdout);
    pcl::io::loadPLYFile( cloud_path, cloud );
    std::cout << "OK, read " << cloud.size() << " points\n";

    std::cout << "minmax..."; fflush(stdout);
    Position min_pt, max_pt;
    pcl::getMinMax3D( cloud, min_pt, max_pt );
    Position div = (max_pt - min_pt);
    div.setConstant( div.head<3>().maxCoeff() );
    std::cout << " min: " << min_pt.transpose() << ", max: " << max_pt.transpose() << ", div: " << div.transpose() << std::endl;

    if ( sceneSize(0) > 0.f )
    {
        div(0) = sceneSize(0) / div(0);
        div(1) = sceneSize(1) / div(1);
        if ( div(2) > 0.f )
            div(2) = sceneSize(2) / div(2);
    }

    Scalar maxDist = 0.;
    if ( doOrigin )
    {
        maxDist = std::max( (max_pt - origin).norm(), (min_pt - origin).norm() );
    }

    float chance = N / (float)cloud.size();
    out_cloud.reserve( N > 0 ? N : cloud.size() );
    for ( size_t pid = 0; pid != cloud.size(); ++pid )
    {
        if ( (N<=0) || ((rand() / (float)RAND_MAX) < chance) )
        {
            auto pnt = cloud.at(pid);

            if ( sceneSize(0) > 0.f )
            {
                pnt.getVector4fMap() -= min_pt;
                pnt.getVector4fMap()(0) *= div(0);
                pnt.getVector4fMap()(1) *= div(1);
                pnt.getVector4fMap()(2) *= div(2);
            }

            if ( doOrigin )
            {
                Eigen::Matrix<int,3,1> colour  = CloudColouring::getColourForDistance( (pnt.getVector3fMap() - origin.head<3>()).norm()
                                                                                     , maxDist, Scalar(0.) );
                pnt.r = colour(0);
                pnt.g = colour(1);
                pnt.b = colour(2);
            }

            out_cloud.push_back( pnt );
        }
    }

    std::cout << "saving..."; fflush(stdout);
    std::stringstream ss;
    ss << cloud_path.substr( 0, cloud_path.find(".ply") ) << "_" << (int)N << ".ply";
    pcl::io::savePLYFile( ss.str(), out_cloud );
    std::cout << " saved to " << ss.str() << std::endl;

    return EXIT_SUCCESS;
}
