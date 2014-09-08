#include "globfit2/visualization/visualization.h"

#include <pcl/console/parse.h>

#include <vector>
#include "globfit2/primitives/linePrimitive2.h"
#include "globfit2/primitives/pointPrimitive.h"

#include "globfit2/io/io.h"                      // readPrimitives()
#include "globfit2/visualization/visualizer.hpp" // GF2::Visualizer
#include "globfit2/util/pcl_util.hpp"            // cloudToVector()

namespace GF2 {
namespace vis {

vis::MyVisPtr
lines::showLines( vis::lines::PrimitiveContainerT  const& lines
                , vis::lines::PointContainerT      const& points
                , vis::lines::Scalar               const  scale
                , Eigen::Matrix<lines::Scalar,3,1> const& colour
                , bool                             const  spin
                , std::vector<lines::Scalar>       const* angles
                , bool                             const  show_ids
                , char                             const  use_tags )
{
    return GF2::Visualizer<vis::lines::PrimitiveContainerT
                          ,vis::lines::PointContainerT>
                          ::show<vis::lines::Scalar>
                          ( lines, points, scale
                          , colour
                          , spin
                          , angles
                          , show_ids
                          , use_tags
                          );
}

int
lines::showLinesCli( int argc, char** argv )
{
    using namespace GF2::vis::lines;
//    typedef GF2::LinePrimitive2                             PrimitiveT;
//    typedef std::vector< std::vector<PrimitiveT> >          PrimitiveContainerT;
//    typedef GF2::PointPrimitive                             PointPrimitiveT;
//    typedef std::vector<PointPrimitiveT>                    PointContainerT;
//    typedef typename PointPrimitiveT::Scalar                Scalar;

    if ( pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
    {
        std::cout << "[" << __func__ << "]: " << "Usage: gurobi_opt --show\n"
                  << "\t--dir \tthe directory containing the files to display\n"
                  << "\t-p,--prims \tthe primitives file name in \"dir\"\n"
                  << "\t--cloud \tthe cloud file name in \"dir\"\n"
                  << "\t[--scale \talgorithm parameter]\n"
                  << "\t[-a,--assoc \tpoint to line associations]\n"
                  << "\t[--no-rel \tdon't show perfect relationships as gray lines]\n"
                  << "\t[--use-tags \tuse associations to create line segments]\n"
                  << "\t[--ids \tshow point GID-s and line GIDs]\n"
                  << "\t[--no-clusters \tdon't show the \"ellipses\"]\n"
                  << "\t[--pop-limit \tpoplation limit for small patches]\n"
                  << "\t[--title \t window title]"
                  << std::endl;
        return EXIT_SUCCESS;
    }

    std::string title = "";
    pcl::console::parse_argument( argc, argv, "--title", title );

    int pop_limit = 10;
    pcl::console::parse_argument( argc, argv, "--pop-limit", pop_limit );

    std::string dir = ".";
    if ( pcl::console::parse_argument( argc, argv, "--dir", dir) < 0 )
    {
        std::cerr << "[" << __func__ << "]: " << "no directory specified by --dir ...assuming local \".\"" << std::endl;
        //return EXIT_FAILURE;
    }
    int err = EXIT_SUCCESS;

    std::string primitives_file;
    if ( (pcl::console::parse_argument( argc, argv, "--prims", primitives_file) < 0) && (pcl::console::parse_argument( argc, argv, "-p", primitives_file) < 0) )
    {
        std::cerr << "[" << __func__ << "]: " << "no primitive file specified by --prims ...exiting" << std::endl;
        return EXIT_FAILURE;
    }

    PrimitiveContainerT lines;
    err = GF2::io::readPrimitives<PrimitiveT, typename PrimitiveContainerT::value_type>( lines, dir + "/" + primitives_file );
    if ( EXIT_SUCCESS != err )
    {
        std::cerr << "[" << __func__ << "]: " << "failed to read " << dir + "/" + primitives_file << "...exiting" << std::endl;
        return err;
    }

    // read cloud
    PointContainerT points;
    std::string cloud_file = "cloud.ply";
    if ( pcl::console::parse_argument( argc, argv, "--cloud", cloud_file) < 0 && !boost::filesystem::exists(cloud_file) )
    {
        std::cerr << "[" << __func__ << "]: " << "no cloud file specified by --cloud ...exiting" << std::endl;
        return EXIT_FAILURE;
    }

    // convert cloud
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud( new pcl::PointCloud<pcl::PointXYZRGB> );
    pcl::io::loadPLYFile( dir + "/" + cloud_file, *cloud );
    GF2::pclutil::cloudToVector<PointPrimitiveT::Allocator>( points, cloud );

    // parse associations

    std::string assoc_file;
    if ( (pcl::console::parse_argument(argc,argv,"--assoc",assoc_file) > 0) || (pcl::console::parse_argument(argc,argv,"-a",assoc_file) > 0) )
    {
        std::vector<std::pair<int,int> > points_primitives;
        std::map<int,int>                linear_indices; // <pid,lid>
        GF2::io::readAssociations( points_primitives, dir + "/" + assoc_file, &linear_indices );

        for ( size_t pid = 0; pid != points_primitives.size(); ++pid )
        {
            if ( points_primitives[pid].first < static_cast<int>(points.size()) )
            {
                points[ pid ].setTag( PointPrimitiveT::GID, points_primitives[pid].first );
            }
            else
                std::cerr << "[" << __func__ << "]: " << "overindexed pid: " << pid << " >= " << points.size() << " points.size(), skipping..." << std::endl;
        }
    }

    // angles
    std::vector<Scalar> angles(3);
    {
        angles[0] = 0;
        angles[1] = M_PI_2;
        angles[2] = M_PI;
    }

    // tags = point colours
    char use_tags = pcl::console::find_switch( argc, argv, "--use-tags" );
    if ( use_tags )
        use_tags += (char)pcl::console::find_switch( argc, argv, "--no-clusters" );

    Scalar scale = 0.1;
    pcl::console::parse_argument(argc, argv, "--scale", scale);

    bool dont_show_rels = pcl::console::find_switch( argc, argv, "--no-rel" );
    bool show_ids       = pcl::console::find_switch( argc, argv, "--ids" );

    GF2::Visualizer<PrimitiveContainerT,PointContainerT>::show<Scalar>( lines
                                                                      , points
                                                                      , scale
                                                                      , (Eigen::Vector3f() << 1,0,0).finished()
                                                                      , /*        spin: */ true
                                                                      , /* connections: */ dont_show_rels ? NULL : &angles
                                                                      , /*    show_ids: */ show_ids
                                                                      , /*    use_tags: */ use_tags
                                                                      , /*   pop-limit: */ pop_limit
                                                                      , /*       title: */ title
                                                                      );
    return EXIT_SUCCESS;
} // ... Solver::show()

} //...ns vis
} //...ns GF2
