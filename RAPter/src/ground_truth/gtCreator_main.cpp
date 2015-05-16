#include "ground_truth/gtCreator.h"
#include "pcl/visualization/pcl_visualizer.h"

#include "pcl/io/ply_io.h"
#include "pcl/console/parse.h"

int main( int argc, char** argv )
{
    using namespace am;
    typedef GF2::GTCreator::PointT PointT;
    srand(123456);

    std::cout << "Hello gtcreator world" << std::endl;
    if ( pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
    {
        std::cout << "usage: " << argv[0] << " --img <image> [--N 200] [--noise 0.005] [--scene_size 1.0]" << std::endl;
        return 0;
    }

    int     n_points = 200;
    pcl::console::parse_argument( argc, argv, "--N", n_points );
    float   noise    = 0.005f;
    pcl::console::parse_argument( argc, argv, "--noise", noise );

    float scene_size = 1.f;
    if ( argc > 4 ) scene_size = atof( argv[4] );

    pcl::PointCloud<PointT>::Ptr cloud( new pcl::PointCloud<PointT>() );
//    if ( argc > 3 )
//        GF2::GTCreator::sampleImage( cloud, std::string(argv[3]), n_points, noise, scene_size );
//    else
//        GF2::GTCreator::run( cloud, "rects", n_points, noise );

    std::cout << "run finished with " << cloud->size() << " points" << std::endl;

    pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
    vptr->setBackgroundColor( .5, .5, .6 );
    vptr->addPointCloud( cloud );
    vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3., "cloud" );

    vptr->spin();

    std::string out_name = "simple.ply";
    if ( argc > 3 )
    {
        std::string             input_str( argv[3]   );
        boost::filesystem::path input    ( input_str );

        std::stringstream ss;
        ss << (input.parent_path() / input.stem()).string() << "_" << n_points << "_" << noise << ".ply";
        out_name = ss.str();
    }

    pcl::io::savePLYFile( out_name, *cloud );
    std::cout << "saved " << out_name << std::endl;

    return EXIT_SUCCESS;
}
