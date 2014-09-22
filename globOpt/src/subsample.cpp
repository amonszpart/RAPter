#include "pcl/console/parse.h"
#include "boost/filesystem.hpp"
#include "pcl/io/io.h"
#include "pcl/io/ply_io.h"
#include "pcl/common/common.h"

int subsample( int argc, char** argv )
{
    typedef Eigen::Vector4f Position;

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

    float scene_size = 1.f;
    if ( pcl::console::parse_argument( argc, argv, "--scene-size", scene_size) < 0 )
    {
        std::cerr << "[" << __func__ << "]: " << "\"--scene-size 1\" assumed to normalize scene to" << std::endl;
    }

    Position sceneSize; sceneSize << scene_size, scene_size, scene_size, 1;

    if ( !valid_input || (pcl::console::find_switch(argc,argv,"-h")) || (pcl::console::find_switch(argc,argv,"--help")) )
    {
        std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << "--subsample\n"
                  << "\t--cloud " << cloud_path << "\n"
                  << "\t--N " << N
                  << "\t--scene-size " << scene_size
                  << "\n";

        return EXIT_FAILURE;
    }

    srand( 123456 );

    pcl::PointCloud<pcl::PointXYZ> cloud, out_cloud;
    out_cloud.reserve( cloud.size() );
    std::cout << "reading " << cloud_path << std::endl;
    pcl::io::loadPLYFile( cloud_path, cloud );

    Position min_pt, max_pt;
    pcl::getMinMax3D( cloud, min_pt, max_pt );
    std::cout << "min: " << min_pt.transpose() << ", max: " << max_pt.transpose() << std::endl;
    Position div = (max_pt - min_pt);
    div.setConstant( div.maxCoeff() );

    if ( scene_size > 0.f )
    {
        div(0) = sceneSize(0) / div(0);
        div(1) = sceneSize(1) / div(1);
        if ( div(2) > 0.f )
            div(2) = sceneSize(2) / div(2);
    }

#if 1
    float chance = N / (float)cloud.size();
    for ( size_t pid = 0; pid != cloud.size(); ++pid )
    {
        if ( (rand() / (float)RAND_MAX) < chance )
        {
            auto pnt = cloud.at(pid);

            if ( scene_size > 0.f )
            {
                pnt.getVector4fMap() -= min_pt;
                pnt.getVector4fMap()(0) *= div(0);
                pnt.getVector4fMap()(1) *= div(1);
                pnt.getVector4fMap()(2) *= div(2);
            }
            std::cout << "adding " << pnt.getVector3fMap().transpose() << " from " << cloud.at(pid).getVector3fMap().transpose() << std::endl;
            out_cloud.push_back( pnt );
        }
    }
#else
    out_cloud.push_back( pcl::PointXYZ(5,5,0) );
    out_cloud.push_back( pcl::PointXYZ(5,6,0) );
    out_cloud.push_back( pcl::PointXYZ(6,5,0) );
    out_cloud.push_back( pcl::PointXYZ(6,6,0) );
#endif

    std::stringstream ss;
    ss << cloud_path.substr( 0, cloud_path.find(".ply") ) << "_" << (int)N << ".ply";
    pcl::io::savePLYFile( ss.str(), out_cloud );

    return EXIT_SUCCESS;
}
