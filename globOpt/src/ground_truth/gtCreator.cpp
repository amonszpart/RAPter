#include "globfit2/ground_truth/gtCreator.h"

#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>

#include "pcl/point_cloud.h"
#include "pcl/point_types.h"
#include "pcl/filters/voxel_grid.h"

#include "globfit2/my_types.h" // pclallocator
#include "globfit2/primitives/linePrimitive.h"
#include "globfit2/ground_truth/ground_truth.h"
#include "globfit2/optimization/candidateGenerator.h"

#if 0
int GF2::GTCreator::run( pcl::PointCloud<PointT>::Ptr &cloud, std::string /*gt_name*/, int gt_nPoints, float gt_noise )
{
    using namespace am;
    using std::vector;

    Eigen::Vector3f sensor_origin( .0f, .0f, .0f );
    const float Z                   = 0.f;
    const float gt_scene_size       = 1.f;
    const bool  gt_clutter_lines    = false;

    // GT
    cv::Mat                       img;
    vector<LinePrimitive>         lines;
    std::vector<float>            angles;

    // code
    {
        //gtImg_stratified( img, lines, angles, 4, 0, 1.f );

        if ( img.empty() )
        {
            img.create( 640, 640, CV_8UC1 );
            img.setTo( 255 );
        }
        std::vector<cv::Point2i> pnts = {
            //{ 87, 125}, {296, 303}
            //, {419,  82}, {597, 243}
            { 39, 443}, {245, 576}
            , {378, 402}, {544, 522}
        };

        for ( size_t i = 0; i != pnts.size(); i+=2 )
            gtRect( img, lines, pnts[i], pnts[i+1],0.f,NULL, 1.f );
    }

    if ( false )
    {
        gtRandomLines( img, lines, gt_clutter_lines, gt_scene_size );
        cv::imshow( "img", img );
    }
    cv::imshow("img",img);
    cv::waitKey( 100 );

    // count GT k
    //desiredK = lines.size();
    cv::imwrite( "input.png", img );

    // IMAGE -> CLOUD
    if ( img.empty() )
    {
        std::cerr << "image input disables...I need gt lines provided..." << std::endl;
        return EXIT_FAILURE;
    }

    const int n_points = gt_nPoints;
    CandidateGenerator::image_2_2DCloud<PCLPointAllocator<3> >( /*     out_cloud: */ *cloud
                                                                , /*        path: */ img
                                                                , /*    N_points: */ n_points
                                                                , /*           Z: */ Z
                                                                , /* scene scale: */ gt_scene_size );

    smartgeometry::addGaussianNoise<MyPoint>( cloud
                                              , NULL
                                              , { gt_noise, gt_noise, 0.f }
                                              , Eigen::Vector3f::Zero()
                                              , &sensor_origin );
    return EXIT_SUCCESS;
}
#endif

int GF2::GTCreator::sampleImage( pcl::PointCloud<PointT>::Ptr &cloud
                                 , std::string                  img_path
                                 , int                          gt_nPoints
                                 , float                        gt_noise
                                 , float                        gt_scene_size
                                 , float                        filter_size
                                 , Eigen::Vector3f              *sensor_origin /* = NULL */ )
{
    const float Z = 0.f;
    //Eigen::Vector3f sensor_origin( .0f, .0f, .0f );
    cv::Mat img = cv::imread( img_path, cv::IMREAD_GRAYSCALE );

    imshow("img", img);
    cv::waitKey( 50 );

    const int n_points = gt_nPoints;
    std::cout << "calling with " << n_points << ", " << gt_scene_size << std::endl;
    CandidateGenerator::image_2_2DCloud<PCLPointAllocator<3> >( /*     out_cloud: */ *cloud
                                                                , /*        path: */ img
                                                                , /*    N_points: */ n_points
                                                                , /*           Z: */ Z
                                                                , /* scene scale: */ gt_scene_size );

    // Create the filtering object
    if ( filter_size > 0.f )
    {
        // convert input
        pcl::PCLPointCloud2::Ptr unfiltered_cloud2( new pcl::PCLPointCloud2() ), filtered_cloud2( new pcl::PCLPointCloud2() );
        pcl::toPCLPointCloud2( *cloud, *unfiltered_cloud2 );

        // filter
        pcl::VoxelGrid<pcl::PCLPointCloud2> sor;
        sor.setInputCloud( unfiltered_cloud2 );
        sor.setLeafSize (filter_size, filter_size, filter_size );
        sor.filter( *filtered_cloud2 );

        pcl::fromPCLPointCloud2( *filtered_cloud2, *cloud );
    }

    smartgeometry::addGaussianNoise<PointT>( cloud
                                           , NULL
                                           , { gt_noise, gt_noise, 0.f }
                                           , Eigen::Vector3f::Zero()
                                           , sensor_origin );
    return EXIT_SUCCESS;
}


