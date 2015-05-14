#include "globfit2/util/parse.h"
#include "globfit2/io/depthIo.hpp"
#include "boost/filesystem.hpp" // exists()
#include "pcl/point_types.h" // pcl::PointXYZRGB
#include "pcl/point_cloud.h" // pcl::PointCloud
#include "pcl/io/ply_io.h" // savePly

int printUsage(int argc, char *argv[])
{
    std::cout << "Usage: " << argv[0] << "\n"
              << "\t--in depth.dat" << "\n"
              << "\t[--rgb colour.png]" << "\n"
              << "\t[-o full_out_path]" << "\n";
    return EXIT_FAILURE;
}

int main(int argc, char **argv )
{
    typedef pcl::PointCloud<pcl::PointXYZRGB> PclCloud;
    int err = EXIT_SUCCESS;

    // read depth
    cv::Mat depth;
    {
        std::string in_path;
        if (    (GF2::console::parse_argument(argc,argv,"--in",in_path) < 0)
             || !boost::filesystem::exists(in_path) )
        {
            return printUsage(argc,argv);
        }
        depth = GF2::io::loadDepth( in_path );
    }

    // read colour
    cv::Mat rgb;
    {
        std::string rgb_path;
        if ( GF2::console::parse_argument(argc,argv,"--rgb",rgb_path) < 0 )
        {
            rgb = cv::imread( rgb_path, cv::IMREAD_UNCHANGED );
        }
    }

    std::string out_path = "./cloud.ply";
    GF2::console::parse_argument( argc,argv,"-o", out_path );

    // convert to cloud
    PclCloud cloud;
    if ( EXIT_SUCCESS == err )
    {
        err = GF2::io::rgbd2PointCloud<ushort>( cloud, depth, cv::Mat(), /* alpha: */ 1/1000.f );
        if ( err != EXIT_SUCCESS )
        {
            std::cerr << "[" << __func__ << "]: " << "rgbd2PointCloud exited with error " << err << std::endl;
        }
    }

    // save cloud
    if ( EXIT_SUCCESS == err )
    {
        pcl::PLYWriter writer;
        writer.write( out_path, cloud, /* binary: */ false, /* use_camera: */ false );
        std::cout << "wrote " << cloud.size() << " points to " << out_path << std::endl;
    }

    return err;
}

