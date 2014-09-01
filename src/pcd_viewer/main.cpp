#include <iostream>
#include <vector>
#include <string>

#include <my_types.h>
#include "pcl/io/pcd_io.h"
#include "pcl/io/ply_io.h"

#include "pcl/visualization/pcl_visualizer.h"

// debug planeprimitive
#include "primitives/planePrimitive.h"
#include "primitives/linePrimitive.h"
#include "lineClustering.hpp"
#include "pcl/point_types.h"

template <typename PrimitiveT> int visualize( int argc, char** argv );

// /home/bontius/Downloads/house_050_noisy00000.pcd /home/bontius/Downloads/house_051_noisy00000.pcd
int main( int argc, char **argv )
{
    using std::vector;
    using std::string;
    using std::cout;
    using std::endl;
    using am::MyCloud;
    using namespace GF2;

    if ( (argc > 2) && (std::string(argv[1]).find("show")!=std::string::npos) )
    {
        //boost::filesystem::path input( std::string(argv[2]) );
        //if ( boost::filesystem::is_directory(input) )

        if ( std::string(argv[2]).find("planes") != std::string::npos )
            return visualize<am::PlanePrimitive>( argc, argv );
        else
            return visualize<am::LinePrimitive>( argc, argv );
    }

    bool project = false;
    if ( (argc > 1) && (std::string(argv[1]).find("project")!=std::string::npos) )
        project = true;

    pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
    vptr->setBackgroundColor( .5, .5, .6 );
    MyCloud::Ptr final_cloud( new MyCloud() );
    for ( int i = project ? 2 : 1; i != argc; ++i )
    {
        MyCloud::Ptr cloud( new MyCloud() );

        if ( !boost::filesystem::exists( argv[i]) ) { std::cerr << argv[i] << " does not exist \n"; continue; }

        pcl::io::loadPCDFile( argv[i], *cloud );

        // clean
        MyCloud::Ptr projected_cloud( new MyCloud() );
        projected_cloud->reserve( cloud->size() );
        for ( size_t pid = 0; pid != cloud->size(); ++pid )
        {
            if ( (cloud->at(pid).x != cloud->at(pid).x)
                 || (cloud->at(pid).y != cloud->at(pid).y)
                 || (cloud->at(pid).z != cloud->at(pid).z) )
                continue;

            projected_cloud->push_back( cloud->at(pid) );
            if ( project )
                projected_cloud->back().z = 0.f;
        }
        cloud = projected_cloud;


        char cloud_name[255];
        sprintf(cloud_name,"cloud%02d",i);
        vptr->addPointCloud( /*(i==1)?tmp_cloud:*/cloud, std::string(cloud_name) );
        vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2., std::string(cloud_name) );
        std::cout << "added " << std::string(argv[i]) << " with " << cloud->size() << " points" << std::endl;

        *final_cloud += *cloud;
    }

    if ( (argc > 1) && final_cloud )
    {
        std::string input_file_str( argv[argc-1] );
        boost::filesystem::path input_file( input_file_str );

        pcl::io::savePLYFile( (input_file.parent_path() / input_file.stem()).string() + "_merged.ply", *final_cloud );
    }

    if ( argc > 1 ) { vptr->spin(); return 0; }

#if 0
    std::cout << "drawing..." << std::endl;
    using am::PlanePrimitive;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud( new pcl::PointCloud<pcl::PointXYZ>() );
    //pcl::PointXYZ pnt(1.f,1.f,1.f);
    //cloud->push_back( pnt );
    vptr->addCoordinateSystem( 0.1, "ref", 0 );

    PlanePrimitive plane( Eigen::Vector3f(1.f,1.f,1.f), Eigen::Vector3f(.1f,.1f,.1f) );
    PlanePrimitive::draw( plane
                          , cloud
                          , pow(3.f,1.f/3.f)
                          , NULL
                          , vptr
                          , std::string("plane0")
                          , 1.,1.,0.
                          , 0);

    Eigen::Vector3f pnt = plane.pos();
    std::cout << "pnt: " << pnt.transpose()
              << ", plane: " << plane().transpose()
              << ", pnt * plane: " << pnt.dot(plane().head<3>()) + plane()(3)
              << std::endl;

    vptr->addArrow( pcl::PointXYZ(pnt(0),pnt(1),pnt(2))
                    , pcl::PointXYZ(0.,0.,0.)
                    , 1., 0., 0., false, "arrow0", 0 );

    vptr->spin();

    cout << "Hello post-opt!" << endl;

    string prims_path = "lines.txt";
    if ( argc > 1 ) prims_path = string( argv[1] );
    string cloud_path = "cloud.pcd";
    if ( argc > 2 ) cloud_path = string( argv[2] );

    // ../../../data/lines.txt
    vector<int> lines_clusters;
    vector<LinePrimitive> prims( am::Primitive<LinePrimitive::Dim>::read<LinePrimitive>(prims_path,&lines_clusters) );
    for ( size_t i = 0; i != prims.size(); ++i )
        cout << "[" << lines_clusters[i] << "] "
                << prims[i]()(0) << ", " << prims[i]()(1) << ", " << prims[i]()(2) << endl;

    MyCloud::Ptr cloud( new MyCloud() );
    pcl::io::loadPCDFile( cloud_path, *cloud );
    for ( size_t pid = 0; pid != cloud->size(); ++pid )
    {
        cout << cloud->at( pid ).x << ", " << cloud->at(pid).y << ", " << cloud->at(pid).z << endl;
    }

    vector<int> points_lines;
    GF2::PostAligner::points2Prims( points_lines
                                    , cloud
                                    , prims
                                    , &(LinePrimitive::point3Distance) );

    int err = GF2::PostAligner::run< vector<LinePrimitive>, am::MyCloud::Ptr, LinePrimitive::Scalar, 2>
            ( prims, cloud, points_lines, lines_clusters, &getLine2DNormal );

    return err;
#endif
    return EXIT_SUCCESS;
}

#include <fstream>

#include "localAnalysis.h"
#include "globFit2.h"

template <class PrimitiveT>
int visualize( int argc, char** argv )
{
    using std::vector;
    //using am::PlanePrimitive;
    typedef vector<PrimitiveT> PrimitivesT;

    if ( !((argc == 5) || (argc == 3))  )
    {
        std::cerr << "[" << __func__ << "]: " << "Usage: " << argv[0] << " --show primitives.txt cloud.pcd solution.txt ...exiting\n";
        return EXIT_FAILURE;
    }

    std::string prims_path(argv[2]), cloud_path, solution_path;
    if ( argc == 3 )
    {
        boost::filesystem::path dir_path = boost::filesystem::path( prims_path ).parent_path();

        boost::filesystem::directory_iterator end_itr; // default construction yields past-the-end
        for ( boost::filesystem::directory_iterator itr( dir_path );
              itr != end_itr;
              ++itr )
        {
            if ( !boost::filesystem::is_directory(itr->status()) )
            {
                if ( itr->path().filename().string().find("cloud") != std::string::npos ) // see below
                    cloud_path = itr->path().string();
                else if ( itr->path().filename().string().find("solution") != std::string::npos ) // see below
                    solution_path = itr->path().string();
            }
        }
    }
    else
    {
        cloud_path = std::string( argv[3] );
        solution_path = std::string( argv[4] );
    }

    // read primitives
    vector<int> lines_clusters;
    PrimitivesT planes = am::Primitive<PrimitiveT::Dim>::template read<PrimitiveT>( prims_path, &lines_clusters );
    std::cout << "[" << __func__ << "]: " << "read " << planes.size() << " planes" << std::endl;

    // read cloud
    am::MyCloud::Ptr cloud( new am::MyCloud() );
    pcl::io::loadPCDFile( cloud_path, *cloud );
    std::cout << "[" << __func__ << "]: " << "read " << cloud->size() << " points" << std::endl;

    // read solution
    am::MaskType opt_mask(planes.size(),0);
    std::vector<float> desired_angles;
    float scale = -1.f;
    {
        std::ifstream f_solution( solution_path );
        if ( !f_solution.is_open() )
        {
            std::cerr << "[" << __func__ << "]: " << "could not open " << solution_path << "...exiting\n";
            return EXIT_FAILURE;
        }

        std::string line;
        enum READ_STATE { PRIMS, SCALE, DESIRED_ANGLES } read_state = PRIMS;
        while ( getline(f_solution, line) )
        {
            // skip comments
            if ( line[0] == '#') continue;

            // read word by word
            std::istringstream  iss( line );
            std::string         word;
            int                 word_id = 0;
            while ( std::getline(iss, word, '\t') )
            {
                // set state
                if ( word_id == 0 )
                {
                    if ( line.find("scale") != std::string::npos )
                        read_state = READ_STATE::SCALE;
                    else if ( line.find("desired_angles") != std::string::npos )
                        read_state = READ_STATE::DESIRED_ANGLES;
                    else
                        read_state = READ_STATE::PRIMS;
                }

                switch ( read_state )
                {
                    case READ_STATE::PRIMS:
                        if ( !word_id ) // first number is the primitive id
                            opt_mask[ atoi(word.c_str()) ] = 1;
                        break;
                    case READ_STATE::SCALE:
                        if ( word_id == 1 )
                            scale = atof(word.c_str());
                        break;
                    case READ_STATE::DESIRED_ANGLES:
                        if ( word_id > 0)
                        {
                            float angle = atof(word.c_str());
                            std::cout << "adding angle " << angle << std::endl; fflush(stdout);
                            desired_angles.push_back( angle );
                            std::cout << "angles is noow ";
                            for(size_t vi=0;vi!=desired_angles.size();++vi)std::cout<<desired_angles[vi]<<" ";std::cout << "\n"; fflush(stdout);
                        }
                        break;
                }
                ++word_id;
            }
        }

        f_solution.close();

        std::cout << "[" << __func__ << "]: " << "read opt_mask: ";
        for ( size_t vi=0; vi != opt_mask.size(); ++vi)
            if ( opt_mask[vi] )
                std::cout << vi << " ";
        std::cout << "\n";
        std::cout << "[" << __func__ << "]: " << "read scale: " << scale << std::endl;
        std::cout << "[" << __func__ << "]: " << "read desired_angles: ";
        for ( size_t vi=0; vi != desired_angles.size(); ++vi)
            std::cout << desired_angles[vi] << " ";
        std::cout << "\n";
    }

    if ( (scale < 0.f) || !desired_angles.size() )
    {
        std::cerr << "[" << __func__ << "]: " << "scale < 0.f or no angles...exiting\n";
        return EXIT_FAILURE;
    }

    auto labels_ptr = (solution_path.find("pearl")!=std::string::npos ? &lines_clusters : NULL );
    //std::cout<<"(*labels_ptr):";for(size_t vi=0;vi!=(*labels_ptr).size();++vi)std::cout<<(*labels_ptr)[vi]<<" ";std::cout << "\\n";

    // hack one line
    {
        //auto it = std::find( opt_mask.begin(), opt_mask.end(), 1 );
        //std::for_each( it+1, opt_mask.end(), []( int &i ) { i = 0; } );
    }

    am::GlobFit2::visualize( planes, opt_mask, cloud, scale, desired_angles, NULL, NULL, labels_ptr, true, false, 5*scale );

    return EXIT_SUCCESS;
}
