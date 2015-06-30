#ifndef __RAPTER_VISUALIZATION_HPP__
#define __RAPTER_VISUALIZATION_HPP__

#include <pcl/console/parse.h>

#include <vector>

#include "rapter/io/io.h"                       // readPrimitives()
#include "rapter/visualization/visualizer.hpp"  // rapter::Visualizer
#include "rapter/util/impl/pclUtil.hpp"         // cloudToVector()
#include "rapter/processing/impl/angleUtil.hpp" // appendAngles...

template <class PrimitiveT>
int
rapter::vis::showCli( int argc, char** argv )
{
    typedef          std::vector< std::vector< PrimitiveT> >              PrimitiveContainerT;
    typedef typename rapter::Visualizer<PrimitiveContainerT,PointContainerT> MyVisT;
    typedef typename MyVisT::Colour                                       Colour;

    if ( pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
    {
        std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << " --show[3D]\n"
                  << "\t--dir \t\t\t The directory containing the files to display\n"
                  << "\t--prims, -p \t\t The primitives file name in \"dir\"\n"
                  << "\t--cloud \t\t The cloud file name in \"dir\"\n"
                  << "\t[ --scale \t\t Algorithm parameter ]\n"
                  << "\t[ --assoc, -a \t\t Point to line associations ]\n\n"

                  << "\t[ --angle-gens \t\t comma separated angle generators. I.e.: 36,90 ]\n"
                  << "\t[ --pop-limit \t\t Poplation limit for small patches ]\n"
                  << "\t[ --point-size \t\t Point render size (6.0) ]\n"
                  << "\t[ --bg-colour \t\t Visualizer background colour i.e. 0.1,0.1,0.1 ]\n"
                  << "\t[ --dir-colours \t colourcode direction IDs' ]\n"
                  << "\t[ --paral-colours \t colourcode parallel directions ]\n"
                  << "\t[ --draw-mode \t\t 0: classic, 1: axis aligned, 2: qhull\n"
                  << "\t              \t\t 4: reprojected, 8: hide primitives, 16: hide unassigned points;  i.e. 28=4+8+16 ]\n"
                  << "\t[ --hull-alpha \t\t QHull alpha parameter ]\n\n"

                  << "\t[ --hide-pts,--no-pts \t Hide point cloud (display only primitives) ]\n"
                  << "\t[ --no-rel \t\t Don't show perfect relationships as gray lines ]\n"
                  << "\t[ --use-tags \t\t Use associations to create line segments ]\n"
                  << "\t[ --ids \t\t Show point GID-s and line GIDs ]\n"
                  << "\t[ --no-clusters \tDon't show the \"ellipses\" ]\n"
                  << "\t[ --pids \t\t Display point ids ]\n"
                  << "\t[ --no-pop \t\t Don't show cluster point counts ]\n"
                  << "\t[ --no-paral \t\t Don't add 0 and 180 to angles - unused ]\n"
                  << "\t[ --no-scale \t\t Hide scale sphere ]\n"
                  << "\t[ --print-angles \t Print angles on relation lines ]\n"
                  << "\t[ --perfect-angle 0.5\t Threshold in degrees, under which a gray line is shown indicating 'perfect relationship' ]\n\n"
                  << "\t[ --show-spatial \t Show spatial distances ]\n"
                  << "\t[ --stretch \t\t Elong primitives by multiplying their extents with this number ]\n"
                  << "\t[ --hide-empty \t\t Don't show empty primitives' ]\n"

                  << "\t[ --title \t\t Window title ]\n"
                  << "\t[ --normals i\t\t show every i-th normal. 0: no normals, 1: every normal ]\n\n"

                  << "\t[ --gids \t\t Show only specific gids (comma separated) ]\n"
                  << "\t[ --dids \t\t Show only specific directions (comma separated) ]\n"
                  << "\t[ --statuses \t\t Show only specific statuses (comma separated) ]\n\n"

                  << "\t[ --save-poly \t\t Save polygons ]\n"
                  << "\t[ --save-hough \t\t Save hough csv]\n"
                  << "\t[ --screenshot \t\t ]\n"
                  << "\t[ --vis-size x,y\t\t ]\n"
                  << std::endl;
        return EXIT_SUCCESS;
    }

    // ------------------

    // show every show_normals-th normal. 0 means no normals, 1 means every normal, 2 for every second normal, etc.
    int show_normals = 0;
    pcl::console::parse_argument( argc, argv, "--normals", show_normals );

    bool    show_pids           = pcl::console::find_switch( argc, argv, "--pids" );
    bool    show_pop            = !pcl::console::find_switch( argc, argv, "--no-pop" );
    bool    print_angles        = pcl::console::find_switch( argc, argv, "--print-angles" );

    float   perfect_angle_limit = 10.e-2;
    pcl::console::parse_argument( argc, argv, "--perfect-angle", perfect_angle_limit );

    bool    paralColours        = pcl::console::find_switch( argc, argv, "--paral-colours" );
    bool    dir_colours         = !paralColours && pcl::console::find_switch( argc, argv, "--dir-colours" );
    bool    hide_points         = pcl::console::find_switch( argc, argv, "--hide-pts" ) || pcl::console::find_switch( argc, argv, "--no-pts" );

    Scalar stretch(1.);
    pcl::console::parse_argument( argc, argv, "--stretch", stretch );

    int     draw_mode           = MyVisT::DRAW_MODE::AXIS_ALIGNED;
    pcl::console::parse_argument( argc, argv, "--draw-mode", draw_mode );

    bool    no_paral            = pcl::console::find_switch( argc, argv, "--no-paral" );
    bool    no_scale_sphere     = pcl::console::find_switch( argc, argv, "--no-scale" );
    bool    save_poly           = pcl::console::find_switch( argc, argv, "--save-poly" );
    bool    save_hough          = pcl::console::find_switch( argc, argv, "--save-hough" );
    bool    show_empty          = !pcl::console::find_switch( argc, argv, "--hide-empty" );


    Scalar  hull_alpha          = 2.;
    pcl::console::parse_argument( argc, argv, "--hull-alpha", hull_alpha);

    std::vector<int> gids;
    pcl::console::parse_x_arguments( argc, argv, "--gids", gids );

    std::vector<int> visSizeVector;
    Eigen::Vector2i visSize;
    if ( pcl::console::parse_x_arguments( argc, argv, "--vis-size", visSizeVector ) >= 0 )
        visSize << visSizeVector.at(0), visSizeVector.at(1);

    std::vector<int> dids;
    pcl::console::parse_x_arguments( argc, argv, "--dids", dids );

    std::vector<int> statuses;
    pcl::console::parse_x_arguments( argc, argv, "--statuses", statuses );

    Colour bg_colour(.1,.1,.1);
    std::vector<float> bg_colours;
    pcl::console::parse_x_arguments( argc, argv, "--bg-colour", bg_colours );
    if ( bg_colours.size() == 3 ) { bg_colour(0) = bg_colours[0]; bg_colour(1) = bg_colours[1]; bg_colour(2) = bg_colours[2]; }

    std::string screenshotPath("");
    pcl::console::parse_argument( argc, argv, "--screenshot", screenshotPath );

    std::string title = "";
    pcl::console::parse_argument( argc, argv, "--title", title );

    bool show_spatial = pcl::console::find_switch( argc, argv, "--show-spatial" );


    int     pop_limit           = 10;
    pcl::console::parse_argument( argc, argv, "--pop-limit", pop_limit );

    std::string dir             = ".";
    if ( pcl::console::parse_argument( argc, argv, "--dir", dir) < 0 )
    {
        std::cerr << "[" << __func__ << "]: " << "no directory specified by --dir ...assuming local \".\"" << std::endl;
    }

    Scalar point_size = 6.;
    pcl::console::parse_argument( argc, argv, "--point-size", point_size );

    // ------------------

    int err = EXIT_SUCCESS;

    std::string primitives_file;
    if ( (pcl::console::parse_argument( argc, argv, "--prims", primitives_file) < 0) && (pcl::console::parse_argument( argc, argv, "-p", primitives_file) < 0) )
    {
        std::cerr << "[" << __func__ << "]: " << "no primitive file specified by --prims ...exiting" << std::endl;
        return EXIT_FAILURE;
    }

    PrimitiveContainerT primitives;
    err = rapter::io::readPrimitives<PrimitiveT, typename PrimitiveContainerT::value_type>( primitives, dir + "/" + primitives_file );
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
    rapter::pclutil::cloudToVector<PointPrimitiveT::Allocator>( points, cloud );
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
        std::vector<std::pair<PidT,LidT> > points_primitives;
        std::map<PidT,LidT>                linear_indices; // <pid,lid>
        rapter::io::readAssociations( points_primitives, dir + "/" + assoc_file, &linear_indices );

        for ( size_t pid = 0; pid != points_primitives.size(); ++pid )
        {
            if ( points_primitives[pid].first < static_cast<int>(points.size()) )
            {
                points[ pid ].setTag( PointPrimitiveT::TAGS::GID, points_primitives[pid].first );
            }
            else
                std::cerr << "[" << __func__ << "]: " << "overindexed pid: " << pid << " >= " << points.size() << " points.size(), skipping..." << std::endl;
        }
    }

    // angles
    AnglesT angle_gens ({Scalar(90.)});
    if ( pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens ) < 0 )
        angle_gens[0] = Scalar(90.);

    AnglesT angles;
    angles::appendAnglesFromGenerators( angles, angle_gens, no_paral, true );

    // tags = point colours
    char use_tags = pcl::console::find_switch( argc, argv, "--use-tags" );
    if ( use_tags )
        use_tags += (char)pcl::console::find_switch( argc, argv, "--no-clusters" );

    Scalar scale = 0.1;
    pcl::console::parse_argument(argc, argv, "--scale", scale);

    bool dont_show_rels = pcl::console::find_switch( argc, argv, "--no-rel" );
    bool show_ids       = pcl::console::find_switch( argc, argv, "--ids" );

    std::set<int> gidsSet;
    {
        gidsSet.insert( gids.begin(), gids.end() );
        if ( gidsSet.size() )
        {
            std::cout<<"gids:";for(size_t vi=0;vi!=gids.size();++vi)std::cout<<gids[vi]<<" ";std::cout << "\n";
        }
    }

    std::set<int>didsSet;
    {
        didsSet.insert( dids.begin(), dids.end() );
        if ( didsSet.size() )
        {
            std::cout<<"dids:";for(size_t vi=0;vi!=dids.size();++vi)std::cout<<dids[vi]<<" ";std::cout << "\n";
        }
    }

    // filter status codes UNSET, ACTIVE and SMALL
    std::set<int> statusSet;
    statusSet.insert( statuses.begin(), statuses.end() );

    // --save-poly will output "cloudRGBNormal" + outStem + ".ply",
    //                     and "plane_mesh"     + outStem + ".ply",
    //                     and "hough"          + outStem + ".ply"
    std::string outStem( "" );
    {
        // Parse ****_itXX.***** to XX in inputPath.
        int iterationFromPath = util::parseIteration( primitives_file );
        // If not, try to find "segments" (regionGrow output), and then "patches" (preMerge output).
        if ( iterationFromPath < 0 )
        {
            if ( primitives_file.find("segments") != std::string::npos )
                outStem = "_segments";
            else if ( primitives_file.find("patches") != std::string::npos )
                outStem = "_patches";
        } //...if iteration number not found
        else
        {
            std::stringstream ss;
            ss << "_it" << iterationFromPath;
            outStem = ss.str();
        }
    } //...outStem parsing

    rapter::Visualizer<PrimitiveContainerT,PointContainerT>::template show<Scalar>( primitives
                                                                               , points
                                                                               , scale
                                                                               , /*               colour: */ (Eigen::Vector3f() << 1,0,0).finished() // probably unused
                                                                               , /*            bg_colour: */ bg_colour
                                                                               , /*                 spin: */ screenshotPath.empty()
                                                                               , /*          connections: */ (dont_show_rels && !show_spatial) ? NULL : &angles
                                                                               , /*             show_ids: */ show_ids
                                                                               , /*             use_tags: */ use_tags
                                                                               , /*            pop-limit: */ pop_limit
                                                                               , /*                title: */ title
                                                                               , /*            show_pids: */ show_pids
                                                                               , /*         show_normals: */ show_normals
                                                                               , /*         show_populat: */ show_pop
                                                                               , /*  perfect_angle_limit: */ perfect_angle_limit
                                                                               , /* print_perfect_angles: */ print_angles
                                                                               , /*          dir_colours: */ dir_colours
                                                                               , /*          hide_points: */ hide_points
                                                                               , /*          filter_gids: */ gidsSet.size  () ? &gidsSet: NULL
                                                                               , /*          filter_dids: */ didsSet.size  () ? &didsSet: NULL
                                                                               , /*        filter_status: */ statusSet.size() ? &statusSet : NULL
                                                                               , /*              stretch: */ stretch
                                                                               , /*            draw_mode: */ draw_mode
                                                                               , /*      no_scale_sphere: */ no_scale_sphere
                                                                               , /*           hull_alpha: */ hull_alpha
                                                                               , /*            save_poly: */ save_poly
                                                                               , /*           point_size: */ point_size
                                                                               , /*           show_empty: */ show_empty
                                                                               , /*         show_spatial: */ show_spatial
                                                                               , /* problemPath (unused): */ ""
                                                                               , /*      parallelColours: */ paralColours
                                                                               , /*              outStem: */ outStem
                                                                               , /*            saveHough: */ save_hough
                                                                               , /*       screenshotPath: */ screenshotPath
                                                                               , /*              visSize: */ visSizeVector.size() ? &visSize : NULL
                                                                               );
    return EXIT_SUCCESS;
} // ... Solver::show()

#endif // __RAPTER_VISUALIZATION_HPP__
