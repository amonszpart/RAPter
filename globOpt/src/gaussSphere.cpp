#include "globfit2/io/inputParser.hpp"
#include "globfit2/util/parse.h"
#include "globfit2/globOpt_types.h"
#include "globfit2/my_types.h"
#include "globfit2/processing/graph.hpp"
#include "globfit2/util/containers.hpp" // class PrimitiveContainer

namespace rapter
{

    /*! \brief Takes oriented, coloured points and saves the normals to a new cloud.
     */
    int gaussSphereCli( int argc, char** argv )
    {
        if ( (argc < 2) || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h")  )
        {
            std::cout << "Usage: " << argv[0] << " cloudRGBNormal.ply [out_name.ply] [--ascii]" << std::endl;
            return EXIT_SUCCESS;
        }
        bool ascii = pcl::console::find_argument( argc, argv, "--ascii" );

        std::string inPath = argv[1], outPath = "";
        if ( !boost::filesystem::exists(inPath) )
        {
            std::cerr << "[" << __func__ << "]: " << "file does not exist, exiting. " << inPath << std::endl;
            return EXIT_FAILURE;
        }

        if ( argc > 2 ) outPath = argv[1];
        else
        {
            boost::filesystem::path inBoostPath( inPath );
            if ( inBoostPath.parent_path().empty() )
                outPath = ".";
            else
                outPath = inBoostPath.parent_path().string();
            outPath += "/" + inBoostPath.stem().string() + "_gaussSphere.ply";
        }
        typedef pcl::PointXYZRGBNormal PclPointT;
        typedef pcl::PointCloud<PclPointT> PclCloudT;
        typedef Eigen::Map< const Eigen::Matrix<float,3,1> > ConstMap3;
        typedef Eigen::Map< Eigen::Matrix<float,3,1> > Map3;

        typename PclCloudT::Ptr cloud( new PclCloudT() );
        pcl::io::loadPLYFile( inPath, *cloud );
        if ( !cloud->size() )
        {
            std::cerr << "[" << __func__ << "]: " << "cloud has zero size, exiting" << std::endl;
            return EXIT_FAILURE;
        }

        typename PclCloudT::Ptr normals( new PclCloudT() );
        normals->reserve( cloud->size() );
        for ( PclCloudT::const_iterator it = cloud->begin(); it != cloud->end(); ++it )
        {
            PclPointT pnt;
            pnt.getVector3fMap() = ConstMap3( it->normal, 3 );
            pnt.rgb              = it->rgb;
            Map3 pntNormal( pnt.normal, 3 );
            pntNormal = it->getVector3fMap();

            normals->push_back( pnt );
        }

        std::cout << "writing to " << outPath << "..."; fflush(stdout);
        pcl::PLYWriter w;
        int err = w.write( outPath, *normals, /* binary_mode: */ !ascii, /* use_camera: */ false );
        std::cout << "...write returned " << err << std::endl;

        return err;
    }
} //...ns GF2

int main( int argc, char *argv[] )
{
    return rapter::gaussSphereCli( argc, argv );

    return 0;
}

