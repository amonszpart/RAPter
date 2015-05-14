#include "pcl/console/parse.h"
#include "boost/filesystem.hpp"
#include "pcl/io/io.h"
#include "pcl/io/ply_io.h"
#include "pcl/common/common.h"
#include "rapter/globOpt_types.h"
#include "rapter/io/io.h"
#include "simpleTypes.h"

int reassign( int argc, char** argv )
{
    typedef Eigen::Vector4f Position;
    typedef rapter::PointContainerT PointContainerT;
    typedef rapter::PointPrimitiveT PointPrimitiveT;
    typedef rapter::_3d::InnerPrimitiveContainerT InnerPrimitiveContainerT;
    typedef rapter::_3d::PrimitiveContainerT PrimitiveContainerT;
    typedef rapter::_3d::PrimitiveT PrimitiveT;

    bool        valid_input             = true;
    std::string cloud_path              = "./cloud.ply", associations_path, input_prims_path;
    float scale = 0.05f;

    // cloud
    if ( (pcl::console::parse_argument( argc, argv, "--cloud", cloud_path) < 0)
         || !boost::filesystem::exists( cloud_path ) )
    {
        std::cerr << "[" << __func__ << "]: " << "--cloud does not exist: " << cloud_path << std::endl;
        valid_input = false;
    }

    // scale
    if ( (pcl::console::parse_argument( argc, argv, "--scale", scale) < 0) && (pcl::console::parse_argument( argc, argv, "-sc", scale) < 0) )
    {
        std::cerr << "[" << __func__ << "]: " << "--scale is compulsory" << std::endl;
        valid_input = false;
    }

    if (    (pcl::console::parse_argument( argc, argv, "-a", associations_path) < 0)
         && (pcl::console::parse_argument( argc, argv, "--assoc", associations_path) < 0)
         && (!boost::filesystem::exists(associations_path)) )
    {
        std::cerr << "[" << __func__ << "]: " << "-a or --assoc is compulsory" << std::endl;
        valid_input = false;
    }

    if (    (pcl::console::parse_argument( argc, argv, "-p", input_prims_path) < 0)
         && (pcl::console::parse_argument( argc, argv, "--prims", input_prims_path) < 0)
         && (!boost::filesystem::exists(input_prims_path)) )
    {
        std::cerr << "[" << __func__ << "]: " << "-p or --prims is compulsory" << std::endl;
        valid_input = false;
    }


    if ( !valid_input )
        return EXIT_FAILURE;

    int err = EXIT_SUCCESS;
    PointContainerT points;
    if ( EXIT_SUCCESS == err )
    {
        err = rapter::io::readPoints<PointPrimitiveT>( points, cloud_path );
        if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
    } //...read points

    typedef std::map<rapter::GidT, InnerPrimitiveContainerT> PrimitiveMapT;
    PrimitiveContainerT planes;
    PrimitiveMapT patches;
    {
        std::cout << "[" << __func__ << "]: " << "reading primitives from " << input_prims_path << "...";
        rapter::io::readPrimitives<PrimitiveT, InnerPrimitiveContainerT>( planes, input_prims_path, &patches );
        std::cout << "reading primitives ok (#: " << planes.size() << ")\n";
    } //...read primitives

    // assign points
    std::cout << "starting assignment" << std::endl; fflush( stdout );
    for ( size_t pid = 0; pid != points.size(); ++pid )
    {
        float min_dist = FLT_MAX, tmp; int min_gid = 0;
        for ( size_t lid = 0; lid != planes.size(); ++lid )
        {
            int gid = -2;
            for ( size_t lid1 = 0; lid1 != planes[lid].size(); ++lid1 )
            {
                if ( !lid1 ) gid = planes[lid][lid1].getTag(PrimitiveT::TAGS::GID);

                tmp = planes[lid][lid1].getDistance( points[pid].pos() );

                if ( (tmp < min_dist) && (tmp < scale) )
                {
                    min_dist = tmp;
                    min_gid = gid;
                }
            }
        }
        //pidGid[ pid ] = min_gid;
        points[pid].setTag( PointPrimitiveT::TAGS::GID, min_gid );
    }
    std::cout << "finishing assignment" << std::endl;
    rapter::io::writeAssociations<PointPrimitiveT>( points, "./points_primitives.schnabel.csv" );




    return EXIT_SUCCESS;
}
