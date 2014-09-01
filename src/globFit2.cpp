#include "globFit2.h"

#include <iostream>
#include <random>
#include <chrono>
#include <omp.h>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "pcl/sample_consensus/sac_model_line.h"
#include "pcl/PointIndices.h"

#include "AMUtil2.h"

#include "kmeans/kMeans.h"
#include "kmeans/anglePrimitive.h"

#include "cuboidRansac.h" // TODO: move to PCLTOOLs

#include "optimization/exhaustiveOptProblem.h"
#include "optimization/simAnnOptProblem.h"

#include "lineClustering.hpp"

// schnabelCandidates
#include "primitives/planePrimitive.h"
#include "schnabelEnv.h"

// MST
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

namespace am
{
    template int
    GlobFit2::optimize( MaskType                                  & min_config
                        , std::vector<LinePrimitive>              const& lines
                        , PrimitiveClustering<std::vector<LinePrimitive> > const* clustering
                        , pcl::PointCloud<MyPoint>::Ptr             cloud
                        , Eigen::Matrix<float,-1,1>          const& lambdas
                        , float                              const  scale
                        , std::vector<float>                 const& desired_angles
                        , int                                const  max_iterations
                        , int                                const  max_step_count
                        , double                             const  threshold         // line inlier threshold
                        , float                              const  trunc_pw_at_angle // in radians
                        , std::vector<int>                   const* indices
                        , int                                const  fixedK
                        , std::vector<int>                        * labels
                        , float                                   * min_e_arg
                        , std::string                        const* p_out_dir );


#warning "DEPR GlobFit2::image_2_2DCloud MOVED TO CANDIDATEGENERATOR"
    /**
     * @brief GlobFit2::image_2_2DCloud Randomly samples a 2D gray image with 3D points at black pixels to create a 3D cloud.
     * @param cloud
     * @param img_path
     * @param N_samples
     * @param Z
     * @return EXIT_SUCCESS/EXIT_FAILURE
     */
    int
    GlobFit2::image_2_2DCloud( pcl::PointCloud<MyPoint>::Ptr &cloud
                               , std::string                 img_path
                               , int                         N_samples
                               , const float                 Z
                               , const float scale )
    {
        return image_2_2DCloud( cloud, cv::imread(img_path,cv::IMREAD_UNCHANGED), N_samples, Z, scale );
    }

    int
    GlobFit2::image_2_2DCloud( pcl::PointCloud<MyPoint>::Ptr  & cloud
                               , cv::Mat                const & img
                               , const int                      N_samples
                               , const float                    Z
                               , const float                    scale )
    {
     //   const float scale = 1.f;

        if ( !cloud )
            cloud = pcl::PointCloud<MyPoint>::Ptr( new pcl::PointCloud<MyPoint>() );

        cv::Mat gray;
        if ( img.channels() > 1 )   cv::cvtColor( img, gray, cv::COLOR_RGB2GRAY );
        else                        gray = img;

        cv::threshold( gray, gray, 200., 255., cv::THRESH_BINARY );
        cv::Mat sampled( gray.clone() ); sampled.setTo(0);

        int   pix_count = gray.cols * gray.rows - cv::countNonZero( gray );
        float prob      = N_samples  / (float)pix_count;
        while ( cloud->size() < static_cast<size_t>(N_samples) )
        {
            pix_count = gray.cols * gray.rows - cv::countNonZero( gray );
            prob      = (N_samples - cloud->size())  / (float)pix_count;

            for ( int y = 0; (y < gray.rows) && (cloud->size() != (size_t)N_samples); ++y  )
            {
                for ( int x = 0; (x < gray.cols) && (cloud->size() != (size_t)N_samples); ++x )
                {
                    float depth = (Eigen::Vector3f() << x, (gray.rows-y), 0).finished().norm();
                    if ( RANDF() > depth / scale / sqrt(2.f) ) continue;
                    if ( (gray.at<uchar>(y,x) < 255) && (RANDF() < prob) )
                    {
                        //std::cout << "cv::sum( gray(cv::Range(y-1,y+1), cv::Range(x-1, x+1))): " << cv::sum( sampled(cv::Range(std::max(0,y-1),std::min(gray.rows,y+1)),
                        //                                                                                             cv::Range(std::max(0,x-1),std::min(gray.cols,x+1))) )[0] << std::endl;
                        if ( cv::sum( sampled(cv::Range(std::max(0,y-1),std::min(gray.rows,y+1)),
                                              cv::Range(std::max(0,x-1),std::min(gray.cols,x+1))) )[0] >= 2 ) continue;

                        MyPoint pnt;
                        pnt.x = x             / (float)gray.cols * scale;
                        pnt.y = (gray.rows-y) / (float)gray.rows * scale;
                        pnt.z = Z;

                        am::util::pcl::setPointColor( pnt, 1, 0, 0 );
                        cloud->push_back( pnt );
                        sampled.at<uchar>(y,x) = 1;
                    }
                }
            }
        }

        cv::imshow( "img", gray );
        cv::waitKey( 20 );

        return EXIT_SUCCESS;
    }


#if 0
    int
    GlobFit2::generateCandidatesSchnabel( std::vector<LinePrimitive>            & lines_out
                                          , pcl::PointCloud<MyPoint>::ConstPtr    cloud
                                          , float                                 Z
                                          , int                                   N_POINT_SAMPLES )
    {
        std::cout << "schnabel..." << std::endl;
        pcl::PointCloud<MyPoint>::Ptr cloud2( new pcl::PointCloud<MyPoint>() );
        const std::vector<float> zs = { /* start: */ Z-.1f, /* stop: */ Z+.1f, /*step:*/ 0.01f };
        for ( size_t pid = 0; pid != cloud->size(); ++pid )
        {
            for ( float z = zs[0]; z < zs[1]; z += zs[2] )
            {
                cloud2->push_back( cloud->at(pid) );
                cloud2->back().z = z;
            }
        }

        std::vector<am::PlanePrimitive> planes;
        SchnabelEnv::run( planes, cloud2, N_POINT_SAMPLES / 10.f );
        if ( !planes.size() )
        {
            std::cerr << "[" << __func__ << "] " << "no PlanePrimitives returned..." << std::endl;
            return EXIT_FAILURE;
        }

        std::vector<am::LinePrimitive> lines; lines.reserve( planes.size() );
        Eigen::Vector4f zPlane; zPlane << 0.f, 0.f, 1.f, 0.f;
        for ( size_t plid = 0; plid != planes.size(); ++plid )
        {
            Eigen::VectorXf line;
            if ( planes[plid].intersectWithPlane( line, zPlane, .01f ) )
                lines.emplace_back( LinePrimitive(line) );
        }

        lines_out.insert( lines_out.end(), lines.rbegin(), lines.rend() );
        std::cout << "schnabel lines inserted: " << lines.size() << std::endl;

        return EXIT_SUCCESS;
    } // generateCandidatesSchnabel()

    /**
     * @brief GlobFit2::showHoughLines Visualizes 2D Hough coordinates.
     * @param vis_params
     * @param lines
     * @param labels
     * @param vptr
     * @param colours
     * @param origin
     * @param viewport
     * @param debug_cloud
     * @return
     */
    int
    GlobFit2::showHoughLines( std::vector<float>                  & vis_params // out: { scale, max_x, max_y }
                              , std::vector<HoughLine>       const& lines
                              , std::vector<int>             const& labels
                              , pcl::visualization::PCLVisualizer::Ptr vptr
                              , std::vector<Eigen::Vector3f> const& colours
                              , Eigen::Vector3f              const& origin
                              , int                                 viewport
                              , pcl::PointCloud<MyPoint>::ConstPtr  debug_cloud )
    {
        float max_r = std::max_element( lines.begin(), lines.end(), [](HoughLine const &a, HoughLine const& b) {
                         return a.R() < b.R();
                     })->R();

        const float scale = 0.5f;
        vis_params = { scale, 2.f*M_PI, max_r };

        pcl::PointCloud<MyPoint>::Ptr cloud( new pcl::PointCloud<MyPoint>() );
        int pid = 0;
        for ( auto const& line : lines )
        {
            MyPoint pnt;
            pnt.x = origin(0) + line.Angle() / vis_params[1] * scale;
            pnt.y = origin(1) + line.R() / vis_params[2] * scale;
            pnt.z = origin(2);
            am::util::pcl::setPointColor( pnt, colours[labels[pid]](0), colours[labels[pid]](1), colours[labels[pid]](2) );

            std::string name = am::util::sprintf( "pid %d", cloud->size() );
            vptr->addText3D( name
                             , pnt
                             , .01
                             , .01,.8,.4
                             , name
                             , 0 );
            cloud->push_back( pnt );

            // debug
            //            if ( debug_cloud )
            //            {
            //                if ( RANDF() < .6 )
            //                {
            //                    vptr->addArrow( pnt
            //                                    , debug_cloud->at(pid)
            //                                    , .5, .6, .9
            //                                    , false,
            //                                    am::util::sprintf("across%d",pid) );
            //                }
            //            }

            ++pid;
        }

        pcl::PointXYZ ang_axis_end(origin(0) + scale, origin(1), origin(2) );
        vptr->addArrow( ang_axis_end
                        , pcl::PointXYZ(origin(0), origin(1), origin(2) )
                        , 0., 0., 0.
                        , false
                        , "ang_axis"
                        , viewport );
        vptr->addText3D( "2PI"
                         , ang_axis_end
                         , .05
                         , 1., 0., 0.
                         , "ang_axis_end_label"
                         , viewport );
        pcl::PointXYZ r_axis_end( origin(0), origin(1) + scale, origin(2) );
        vptr->addArrow( r_axis_end
                        , pcl::PointXYZ(origin(0), origin(1), origin(2) )
                        , 0., 0., 0.
                        , false
                        , "r_axis"
                        , viewport );
        vptr->addText3D( am::util::sprintf("%6.3f", max_r)
                         , r_axis_end
                         , .05
                         , 1., 0., 0.
                         , "r_axis_end_label"
                         , viewport );
        vptr->addCoordinateSystem( .1, "coord_system", viewport );
        vptr->addPointCloud( cloud, "hough_lines", viewport );
        vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3., "hough_lines" );

        return EXIT_SUCCESS;
    }

    /**
     * @brief GlobFit2::cluster Calls kmeans to cluster hough lines with fixed K sizes
     * @param primitives
     * @param labels
     * @param hough_lines
     * @param Ks
     * @param vptr
     * @param debug_cloud
     * @return
     */
    int
    GlobFit2::cluster(
            typename kmeans::Primitive<std::vector<HoughLine>,HoughLine>::Vec & primitives
            , std::vector<int>                                                & labels
            , std::vector<HoughLine>                                     const& hough_lines
            , std::vector<int>                                           const& Ks
            , pcl::visualization::PCLVisualizer::Ptr                            vptr
            , pcl::PointCloud<MyPoint>::ConstPtr                                debug_cloud
            )
    {
        labels.resize( hough_lines.size(), 0 );

        // cluster ANGLE
        kmeans::Primitive<std::vector<HoughLine>,HoughLine>::Vec angle_primitives;
        std::vector<int>                                 angle_labels;
        {
            // init primitives
            std::vector<int> angle_dims = { 1 };
            for ( int k = 0; k != Ks[0]; ++k )
                angle_primitives.emplace_back( new kmeans::AnglePrimitive<std::vector<HoughLine>,HoughLine>(angle_dims) ); // delete in primitives

            // run
            kmeans::KMeans<std::vector<HoughLine>,HoughLine>::cluster( /*   [io]     clusters: */ angle_primitives
                                                                       , /* [io]       labels: */ angle_labels
                                                                       , /*            points: */ hough_lines
                                                                       , /*        dim_coeffs: */ NULL
                                                                       , /*               kpp: */ true );
            // copy primitives to clusters
            for ( size_t k = 0; k != angle_primitives.size(); ++k )
            {
                std::cout << "cluster[" << k << "] size: " << angle_primitives[k]->size() << std::endl;
            }
        }

        // JSON
        std::vector<Eigen::Vector3f> colours = am::util::nColoursEigen( angle_labels.size() );
        kmeans::KMeans<std::vector<HoughLine>,HoughLine>::toJson( "angles.json"
                                                                  , hough_lines.data(), hough_lines.size()
                                                                  , angle_labels
                                                                  , angle_primitives
                                                                  , colours );
        // cluster DIST
        kmeans::Primitive<std::vector<HoughLine>,HoughLine>::Vec tmp_primitives;
        std::vector<int>                                 tmp_labels;
        for ( int prim_id = 0; (size_t)prim_id != angle_primitives.size(); ++prim_id )
        {
            // select points from angle_cluster
            std::vector<HoughLine> tmp_lines;
            for ( int pid = 0; pid != (int)hough_lines.size(); ++pid )
            {
                if ( angle_labels[pid] == prim_id )
                    tmp_lines.push_back( hough_lines[pid] );
            }

            // prepare primitives
            SAFE_CLEAR_VECTOR( tmp_primitives );
            for ( int k = 0; k != Ks[1]; ++k )
                tmp_primitives.emplace_back( new kmeans::PointPrimitive<std::vector<HoughLine>,HoughLine>() ); // delete in primitives

            // select first dimension
            HoughLine::VectorType dim_coeffs( HoughLine::Dim, 1 );
            dim_coeffs << 1.f, 0.f;

            // work
            kmeans::KMeans<std::vector<HoughLine>,HoughLine>::cluster( tmp_primitives
                                                                       , tmp_labels
                                                                       , tmp_lines
                                                                       , &dim_coeffs
                                                                       , /* kpp: */ true );

            // write back labels
            int tmp_line_id = 0;
            for ( int pid = 0; pid != (int)hough_lines.size(); ++pid )
            {
                if ( angle_labels[pid] == prim_id )
                    labels[pid] = primitives.size() + tmp_labels[tmp_line_id++];
            }

            // JSON
            kmeans::KMeans<std::vector<HoughLine>,HoughLine>::toJson( am::util::sprintf("points%d.json",prim_id)
                                                                      , tmp_lines.data(), tmp_lines.size()
                                                                      , tmp_labels
                                                                      , tmp_primitives
                                                                      , colours );

            // copy out
            for ( size_t k = 0; k != tmp_primitives.size(); ++k )
            {
                std::cout << "cluster[" << k << "] size: " << tmp_primitives[k]->size() << std::endl;

                // output
                primitives.push_back( tmp_primitives[k] ); // ownership now in primitives
            }
            tmp_primitives.clear(); // ownership in primitives
        }

        // export
        // JSON
        for ( size_t i = 0; i != angle_primitives.size(); ++i )
        {
            primitives.push_back( angle_primitives[i] );
            primitives.back()->asPointPrimitive()->Pnt()(0) = 0.f; // on x axis
        }
        angle_primitives.clear(); // ownership in primitives

        kmeans::KMeans<std::vector<HoughLine>,HoughLine>::toJson( "points.json"
                                                                  , hough_lines.data(), hough_lines.size()
                                                                  , labels
                                                                  , primitives
                                                                  , colours );

        if ( primitives.back() ) { delete primitives.back(); primitives.back() = NULL; } primitives.pop_back();
        if ( primitives.back() ) { delete primitives.back(); primitives.back() = NULL; } primitives.pop_back();

        // clenaup
        // ownersip at caller now...
        // SAFE_CLEAR_VECTOR( primitives );

        return EXIT_SUCCESS;
    }

    void mouseEventOccurred (const pcl::visualization::MouseEvent &event,
                             void* viewer_void)
    {
        //std::cout << "event.getKeybooard:" << event.getKeyboardModifiers() << std::endl;

        boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = *static_cast<boost::shared_ptr<pcl::visualization::PCLVisualizer> *> (viewer_void);
        if ( event.getButton () == pcl::visualization::MouseEvent::LeftButton
             && event.getType () == pcl::visualization::MouseEvent::MouseButtonRelease
             && !(event.getKeyboardModifiers() & pcl::visualization::KeyboardEvent::Shift) )
        {
            std::cout << "Left mouse button released at position (" << event.getX () << ", " << event.getY () << ")" << std::endl;

            char str[512];
            static int text_id = 0;
            sprintf (str, "text#%03d", text_id++ );
            viewer->addText ("clicked here", event.getX (), event.getY (), str);
        }
    }

    struct VisualizationState
    {
            VisualizationState()
                : clicked_points_3d( new pcl::PointCloud<MyPoint>() )
            {
            }

            // structure used to pass arguments to the callback function
            pcl::PointCloud<MyPoint>::Ptr          clicked_points_3d;
            pcl::visualization::PCLVisualizer::Ptr viewerPtr;
            std::vector<pcl::ModelCoefficients>    circles;
    };

    void
    pp_callback (const pcl::visualization::PointPickingEvent& event, void* args)
    {

        struct VisualizationState* data = (struct VisualizationState *)args;
        if (event.getPointIndex () == -1)
            return;

        MyPoint current_point;
        event.getPoint(current_point.x, current_point.y, current_point.z);
        data->clicked_points_3d->clear();
        data->clicked_points_3d->points.push_back(current_point);

        data->viewerPtr->removeShape( "clicked_name" );
        data->viewerPtr->addText3D( am::util::sprintf( "pid %d", event.getPointIndex() )
                                    , current_point
                                    , .01
                                    , 0., 0., 1.
                                    , "clicked_name"
                                    , 0 );

        // todo: change circles to hough_lines
        //int cid = 0;
        float min_diff = HUGE_VAL;
        for ( size_t i = 1; i != data->circles.size(); ++i )
        {
            Eigen::Vector2f c; c << data->circles[i].values[0], data->circles[i].values[1];
            auto diff = ( c - current_point.getVector3fMap().head<2>() ).norm();
            if ( diff < min_diff )
            {
                min_diff = diff;
                //cid = i;
            }
        }

        // Draw clicked points in red:
        pcl::visualization::PointCloudColorHandlerCustom<MyPoint> red (data->clicked_points_3d, 255, 0, 0);
        data->viewerPtr->removePointCloud("clicked_points");
        data->viewerPtr->addPointCloud(data->clicked_points_3d, red, "clicked_points");
        data->viewerPtr->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "clicked_points");
        std::cout << current_point.x << " " << current_point.y << " " << current_point.z << std::endl;
    }

    int
    GlobFit2::generateCandidates( std::vector<LinePrimitive>           & lines
                                  , pcl::PointCloud<MyPoint>::ConstPtr   cloud
                                  , const float                          threshold )
    {
        std::cerr << "[" << __func__ << "] " << ": deprecated!" << std::endl;
#if 0
        // display
        VisualizationState vis_state;
        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
        {
            vptr->setBackgroundColor( .4, .8, .4 );
            vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3. );
            //vptr->registerMouseCallback( mouseEventOccurred, reinterpret_cast<void*>(&vptr) );
            //vptr->registerPointPickingCallback( &pointPickingCallback, reinterpret_cast<void*>(&vptr) );
            vptr->registerPointPickingCallback( pp_callback, (void*)&vis_state );
            vis_state.viewerPtr = vptr;
        }

        vptr->addPointCloud( cloud );

        // get local neighs
        std::vector<pcl::PointIndices::Ptr> neighbour_indices;
        pcl::PointIndices::Ptr              indices           = am::util::pcl::allIndicesOf( cloud );
        {
            CuboidRansac::getNeighbourhoodIndices( neighbour_indices
                                                   , cloud
                                                   , indices );
        }

        // get line models
        std::vector<pcl::ModelCoefficients::Ptr> models;
        {
            CuboidRansac::getLocalModels( /* [out]         models: */ models
                                          , /*      neighbourhood: */ neighbour_indices
                                          , /*          cloud_ptr: */ cloud
                                          , /*         model type: */ pcl::SACMODEL_LINE
                                          , /*         distThresh: */ threshold );
        }

        // line normals
        //Eigen::Vector3f plane_normal( 0, 0, 1.f ); // perp to xy plane
        //MY_ASSERT( models.size() == indices->indices.size() );
        //for ( size_t i = 0; i < models.size(); ++i )
        //{
        //Eigen::Vector3f line_dir( models[i]->values[3], models[i]->values[4], models[i]->values[5] );
        //MyPoint         pnt = cloud->at( indices->indices[i] );
        //Eigen::Vector3f normal = -1.f*line_dir.cross(plane_normal).normalized();
        //std::string     arrow_name = am::util::sprintf( "line_normal%d", i );
        //            vptr->addArrow( am::util::pcl::asPointXYZ(pnt.getVector3fMap() + normal * 0.1f)
        //                            , pnt
        //                            , .7, .1, .1
        //                            , false
        //                            , arrow_name
        //                            , 0 );
        //}

        // HoughLines - to parameter space
        std::vector<HoughLine> hough_lines( models.size() );
        {
            for ( size_t i = 0; i != models.size(); ++i )
                hough_lines[i] = HoughLine::from3DSpace( models[i]->values );
        }

        // Cluster
        kmeans::Primitive<std::vector<HoughLine>,HoughLine>::Vec primitives;
        std::vector<int>                                 labels;
        cluster( primitives
                 , labels
                 , hough_lines
                 , { 4, 4 }
                 , vptr
                 , cloud );

        //std::vector<Line> lines;
        lines.clear(); lines.reserve(primitives.size());
        for ( size_t k = 0; k < primitives.size(); ++k )
        {
            lines.emplace_back( LinePrimitive(HoughLine::to3DSpace<float>(primitives[k]->asPointPrimitive()->Pnt())) );
        }

        // draw line candidates
        {
            for ( size_t k = 0; k < lines.size(); ++k )
            {
                //                pcl::ModelCoefficients line_coeffs;
                //                line_coeffs.values = lines[k].coeffsVector();

                //                pcl::PointIndices::Ptr        inliers( new pcl::PointIndices() );
                //                pcl::PointCloud<MyPoint>::Ptr on_line_cloud( new pcl::PointCloud<MyPoint>() );
                //                {
                //                    pcl::SampleConsensusModelLine<MyPoint> sacline( cloud );
                //                    if ( indices )
                //                        sacline.setIndices( indices->indices );
                //                    sacline.selectWithinDistance( lines[k].coeffs(), threshold        , inliers->indices );
                //                    sacline.projectPoints       ( inliers->indices , lines[k].coeffs(), *on_line_cloud, false );
                //                }

                //                // select min and max point
                //                if ( !on_line_cloud->size() )
                //                    continue;

                //                float min_dist = 0.f, max_dist = 0.f;
                //                int   min_id   = 0  , max_id = 0;
                //                Eigen::Vector3f p0       = on_line_cloud->at(0).getVector3fMap();
                //                Eigen::Vector3f line_dir = lines[k].dir(); line_dir.normalize();
                //                for ( size_t point_id = 1; point_id != on_line_cloud->size(); ++point_id )
                //                {
                //                    Eigen::Vector3f p1 = on_line_cloud->at(point_id).getVector3fMap();
                //                    Eigen::Vector3f p0p1 = p1-p0;
                //                    float dist = p0p1.dot( p0 + line_dir );
                //                    if ( dist < min_dist )
                //                    {
                //                        min_dist = dist;
                //                        min_id = point_id;
                //                    }
                //                    else if ( dist > max_dist )
                //                    {
                //                        max_dist = dist;
                //                        max_id = point_id;
                //                    }
                //                }

                //vptr->addLine( line_coeffs, am::util::sprintf("hough_line2_%d",k), 0 );

                std::vector<MyPoint> minMax;
                int err = lines[k].getExtent<MyPoint>( minMax
                                                       , cloud
                                                       , threshold
                                                       , indices?&(indices->indices):NULL );
                if ( err != EXIT_SUCCESS ) continue;

                vptr->addLine( am::util::pcl::asPointXYZ(minMax[0].getVector3fMap()
                               - (minMax[1].getVector3fMap()
                               - minMax[0].getVector3fMap())
                        * .5f),
                        am::util::pcl::asPointXYZ(minMax[1].getVector3fMap()
                        - (minMax[0].getVector3fMap()
                        - minMax[1].getVector3fMap())
                        * .5f),
                        0., .3, 1.,
                        am::util::sprintf("hough_line%d",k), 0 );
                //vptr->spin();
                //                vptr->addArrow( on_line_cloud->at(0)
                //                                , pcl::PointXYZ(circles[k].values[0],circles[k].values[1],0.f)
                //                                , .6, .6, .9
                //                                , false
                //                                , am::util::sprintf("centroid2line%d",k)
                //                                , 0 );
            }
        }

        // vis_params
        std::vector<float>    vis_params; // { scale, x_normalizer, y_normalizer }
        const Eigen::Vector3f origin( 1.f, 0.f, 0.f );
        std::vector<Eigen::Vector3f> colours = am::util::nColoursEigen( primitives.size() );

        // draw points
        showHoughLines( vis_params
                        , hough_lines
                        , labels
                        , vptr
                        , colours
                        , origin
                        , 0
                        , /* debug: */ cloud );

        // draw circles
        vis_state.circles.resize( primitives.size() );
        {
            for ( size_t k = 0; k < vis_state.circles.size(); ++k )
            {
                float r = 0.f;
                for ( size_t pid = 0; pid != hough_lines.size(); ++pid )
                {
                    if ( labels[pid] == static_cast<int>(k) )
                    {
                        Eigen::Vector2f centroid = primitives[k]->asPointPrimitive()->Pnt();
                        Eigen::Vector2f dist( Eigen::Vector2f::Zero() );
                        // r dist:
                        dist(0) = (centroid(0) - hough_lines[pid].pos()(0)) / vis_params[2] * vis_params[0];
                        // ang dist:
                        dist(1) = (centroid(1) - hough_lines[pid].pos()(1)) / vis_params[1] * vis_params[0];
                        r = std::max( r, dist.norm() );
                    }
                }
                vis_state.circles[k].values = { origin(0) + primitives[k]->asPointPrimitive()->Pnt()(1) / vis_params[1] * vis_params[0],
                                                origin(1) + primitives[k]->asPointPrimitive()->Pnt()(0) / vis_params[2] * vis_params[0],
                                                r };
                vptr->addCircle( vis_state.circles[k], am::util::sprintf("circle%d",k),0 );
            }
        }

        vptr->spin();

        // ownership at caller
        SAFE_CLEAR_VECTOR( primitives );
#endif
        return EXIT_SUCCESS;
    }
#endif

} // ns am
