#ifndef GF2_VISUALIZATION_HPP
#define GF2_VISUALIZATION_HPP

#include <pcl/console/parse.h>

#include <vector>

#include "globfit2/io/io.h"                      // readPrimitives()
#include "globfit2/visualization/visualizer.hpp" // GF2::Visualizer
#include "globfit2/util/pcl_util.hpp"            // cloudToVector()

template <class PrimitiveT>
int
GF2::vis::showCli( int argc, char** argv )
{
    typedef std::vector< std::vector< PrimitiveT> > PrimitiveContainerT;

    if ( pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
    {
        std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << " --show[3d]\n"
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
                  << "\t[--pids \tdisplay point ids]\n"
                  << "\t[--title \t window title]\n"
                  << "\t[--normals i\t show every i-th normal. 0: no normals, 1: every normal]\n"
                  << "\t[--angle-gens \t comma separated angle generators. I.e.: 36,90]\n"
                  << "\t[--no-pop \t Don't show cluster point counts]\n"
                  << "\t[--print-angles \t Print angles on relation lines]\n"
                  << "\t[--perfect-angle 0.5\t Threshold in degrees, under which a gray line is shown indicating 'perfect relationship']\n"
                  << std::endl;
        return EXIT_SUCCESS;
    }

    // show every show_normals-th normal. 0 means no normals, 1 means every normal, 2 for every second normal, etc.
    int show_normals = 0;
    pcl::console::parse_argument( argc, argv, "--normals", show_normals );

    bool show_pids = pcl::console::find_switch( argc, argv, "--pids" );
    bool show_pop = !pcl::console::find_switch( argc, argv, "--no-pop" );
    bool print_angles = pcl::console::find_switch( argc, argv, "--print-angles" );
    float perfect_angle_limit = 10.e-5;
    pcl::console::parse_argument( argc, argv, "--perfect-angle", perfect_angle_limit );

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
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud( new pcl::PointCloud<pcl::PointNormal> );
    pcl::io::loadPLYFile( dir + "/" + cloud_file, *cloud );
    GF2::pclutil::cloudToVector<PointPrimitiveT::Allocator>( points, cloud );
    for ( int pid = 0; pid != cloud->size(); ++pid )
    {
        points[pid].coeffs()(3) = cloud->at(pid).normal_x;
        points[pid].coeffs()(4) = cloud->at(pid).normal_y;
        points[pid].coeffs()(5) = cloud->at(pid).normal_z;
    }

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
    std::vector<Scalar> angle_gens(1);
    if ( pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens ) < 0 )
        angle_gens[0] = Scalar(90.);

    std::vector<Scalar> angles;
    processing::appendAnglesFromGenerators( angles, angle_gens, true );
//    {
//        angles[0] = 0;
//        angles[1] = M_PI_2;
//        angles[2] = M_PI;
//    }

    // tags = point colours
    char use_tags = pcl::console::find_switch( argc, argv, "--use-tags" );
    if ( use_tags )
        use_tags += (char)pcl::console::find_switch( argc, argv, "--no-clusters" );

    Scalar scale = 0.1;
    pcl::console::parse_argument(argc, argv, "--scale", scale);

    bool dont_show_rels = pcl::console::find_switch( argc, argv, "--no-rel" );
    bool show_ids       = pcl::console::find_switch( argc, argv, "--ids" );

    GF2::Visualizer<PrimitiveContainerT,PointContainerT>::template show<Scalar>( lines
                                                                               , points
                                                                               , scale
                                                                               , /*               colour: */ (Eigen::Vector3f() << 1,0,0).finished()
                                                                               , /*                 spin: */ true
                                                                               , /*          connections: */ dont_show_rels ? NULL : &angles
                                                                               , /*             show_ids: */ show_ids
                                                                               , /*             use_tags: */ use_tags
                                                                               , /*            pop-limit: */ pop_limit
                                                                               , /*                title: */ title
                                                                               , /*            show_pids: */ show_pids
                                                                               , /*         show_normals: */ show_normals
                                                                               , /*         show_populat: */ show_pop
                                                                               , /*  perfect_angle_limit: */ perfect_angle_limit
                                                                               , /* print_perfect_angles: */ print_angles
                                                                               );
    return EXIT_SUCCESS;
} // ... Solver::show()

#endif // GF2_VISUALIZATION_HPP
