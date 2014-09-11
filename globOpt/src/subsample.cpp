#include "pcl/console/parse.h"
#include "boost/filesystem.hpp"
#include "pcl/io/io.h"
#include "pcl/io/ply_io.h"

int subsample( int argc, char** argv )
{
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

    if ( !valid_input || (pcl::console::find_switch(argc,argv,"-h")) || (pcl::console::find_switch(argc,argv,"--help")) )
    {
        std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << "--subsample\n"
                  << "\t--cloud " << cloud_path << "\n"
                  << "\t--N " << N << "\n";

        return EXIT_FAILURE;
    }

    srand( 123456 );

    pcl::PointCloud<pcl::PointXYZ> cloud, out_cloud;
    out_cloud.reserve( cloud.size() );
    std::cout << "reading " << cloud_path << std::endl;
    pcl::io::loadPLYFile( cloud_path, cloud );
    float chance = N / (float)cloud.size();
    for ( size_t pid = 0; pid != cloud.size(); ++pid )
    {
        if ( (rand() / (float)RAND_MAX) < chance )
        {
            out_cloud.push_back( cloud.at(pid) );
        }
    }

    std::stringstream ss;
    ss << cloud_path.substr( 0, cloud_path.find(".ply") ) << "_" << (int)N << ".ply";
    pcl::io::savePLYFile( ss.str(), out_cloud );

    return EXIT_SUCCESS;
}
