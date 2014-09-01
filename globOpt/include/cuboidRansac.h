#ifndef CUBOIDRANSAC_H
#define CUBOIDRANSAC_H

#include <string>
#include "my_types.h"
#include "pcl/point_cloud.h"
#include "pcl/PointIndices.h"
#include "pcl/ModelCoefficients.h"
#include "pcl/sample_consensus/model_types.h"

#include "pcltools/util.hpp" // smartgeometry

namespace am
{

    // this is a really bad name...sorry
    class CuboidRansac
    {
        public:
            /**
             * @brief CuboidRansac::run Takes a pointcloud (from disk), finds a plane, projects all points to that plane, and finds a 2D concave hull,
             *                          selects points close to that hull, and finds local line suggestions based on these points
             * @param cloud_path        Source for the 3D pointcloud.
             * @param alpha             How smooth you want the pointcloud (I think the lower, the less smooth...but don't know)
             * @param MyPointT          Any pcl pointtype. (pcl::PointXYZRGB)
             * @return EXIT_SUCCESS
             */
            template <typename MyPointT>
            static int
            run( boost::shared_ptr<pcl::PointCloud<MyPointT> >   &hull_cloud
                 , std::string                                    cloud_path
                 , double                                         alpha      = .05 );

            // I moved this to PCL tools since...but anyway, this will return a list of list of ints as the closest K points
            // in the new version, there's a radius limit as well, cause here, it's not guaranteed to be close...
            static int
            getNeighbourhoodIndices(std::vector<pcl::PointIndices::Ptr>  &neighbour_indices
                                    , pcl::PointCloud<MyPoint>::ConstPtr  cloud
                                    , pcl::PointIndices::ConstPtr         indices   = pcl::PointIndices::ConstPtr()
                                    , std::vector<std::vector<float> >   *distances = NULL
                                    , int                                 K         = 15                            );

            /**
             * @brief getLocalModels          Gives you ransac lines from neighbour indices
             * @param [out] models            List of ransaced lines, one for every entry in neighbour_indices/cloud
             * @param [in ] neighbour_indices List of list of neighbour indices, one list for each point in the cloud.
             * @param cloud                   3D pointcloud
             * @param sac_model               Desired sacmodel, right now only SACMODEL_LINE
             * @param distanceThreshold       Inlier distance threshold for the SAC
             * @return                        EXIT_SUCCESS
             */
            static int
            getLocalModels( std::vector<pcl::ModelCoefficients::Ptr>    & models
                           , std::vector<std::vector<int> >       const & neighbour_indices
                           , pcl::PointCloud<MyPoint>::ConstPtr           cloud
                           , int                                          sac_model         = pcl::SACMODEL_LINE
                           , double                                       distanceThreshold = 0.01                          );
    };

} // ns am

#ifndef __CUBOIDRANSAC_HPP__
#define __CUBOIDRANSAC_HPP__

#include "pcl/point_types.h"
#include "pcl/PointIndices.h"
#include "pcl/io/pcd_io.h"
#include "pcl/filters/project_inliers.h"
#include "pcl/segmentation/sac_segmentation.h"
#include "pcl/PointIndices.h"
#include "pcl/search/kdtree.h"
#include "pcl/surface/concave_hull.h"
#include "pcl/visualization/pcl_visualizer.h"

#include "pcltools/primitives.h" // smartgeom


namespace am
{
    inline Eigen::VectorXf
    coeffs2Vector( pcl::ModelCoefficients::Ptr coeffs )
    {
        Eigen::VectorXf out( coeffs->values.size() );
        std::copy( coeffs->values.begin(), coeffs->values.end(), out.data() );

        return out;
    }

    template <typename MyPointT>
    int
    CuboidRansac::run( boost::shared_ptr<pcl::PointCloud<MyPointT> >   & hull_cloud
                       , std::string                                     cloud_path
                       , double                                          alpha     )
    {
        // read
        typename pcl::PointCloud<MyPointT>::Ptr cloud_ptr ( new pcl::PointCloud<MyPointT>() );
        {
            pcl::io::loadPCDFile( cloud_path, *cloud_ptr );
        }

        // prepare output
        if ( !hull_cloud )  hull_cloud.reset( new pcl::PointCloud<MyPointT>() );

        // select all indices
        //pcl::PointIndices::Ptr indices_ptr( am::util::pcl::allIndicesOf(cloud_ptr) );

        // get plane
        pcl::ModelCoefficients::Ptr plane_coeffs   ( new pcl::ModelCoefficients() );
        pcl::PointIndices::Ptr      inliers_indices( new pcl::PointIndices()      );
        smartgeometry::Primitives::extractPlane<MyPointT>( plane_coeffs, cloud_ptr, &inliers_indices, nullptr );

        // project to plane
        typename pcl::PointCloud<MyPointT>::Ptr cloud_projected( new pcl::PointCloud<MyPointT>() );
        {
            pcl::ProjectInliers<MyPointT> proj;
            proj.setModelType         ( pcl::SACMODEL_PLANE );
            proj.setInputCloud        ( cloud_ptr           );
            proj.setModelCoefficients ( plane_coeffs        );
            proj.filter               ( *cloud_projected    );
        }

        // init visualization
        pcl::visualization::PCLVisualizer::Ptr vptr(  new pcl::visualization::PCLVisualizer() );
        int vp0, vp1;
        vptr->setSize( 1400, 900 );
        vptr->createViewPort( 0, 0, 0.5, 1, vp0 );
        vptr->setBackgroundColor( .8, .8, .8, vp0 );
        vptr->createViewPort( 0.5, 0, 1, 1, vp1 );
        vptr->setBackgroundColor( .8, .8, .8, vp1 );
        vptr->addPointCloud( cloud_projected, "cloud", vp0 );
        vptr->addPlane( *plane_coeffs, 0,0,0,"plane", vp0 );

        // qhull - takes the
        {
            // Create a Convex Hull representation of the projected inliers
            pcl::ConcaveHull<MyPointT>              concave_hull;                                   // object
            typename pcl::PointCloud<MyPointT>::Ptr cloud_hull( new pcl::PointCloud<MyPointT>() );  // output points
            std::vector<pcl::Vertices>              polygons;                                       // output list indexing the points from cloud_hull, in 2D this is size 1
            concave_hull.setAlpha( alpha );                                                         // smoothness
            concave_hull.setInputCloud( cloud_projected );
            concave_hull.reconstruct  ( *cloud_hull, polygons );                                    // work

            // debug output
            std::cerr << "Convex hull has: " << cloud_hull->points.size () << " data points "
                      << " in " << polygons.size() << " polygons " << std::endl;

            // for every border polygon
            for ( size_t poly_id = 0; poly_id != polygons.size(); ++poly_id )
            {
                // for every vertex in polygon boundary - create a line segment using the prev point in the hull,
                // and get points from the cloud close to this segment
                pcl::PointCloud<pcl::PointXYZ>::Ptr polygon_cloud_ptr( new pcl::PointCloud<pcl::PointXYZ>() );
                for ( size_t vert_id = 0; vert_id != polygons[poly_id].vertices.size(); ++vert_id )
                {
                    // line endpoint
                    MyPointT pnt1 = cloud_hull->at( polygons[poly_id].vertices[vert_id] );
                    // copy to polygon_cloud
                    polygon_cloud_ptr->push_back( pcl::PointXYZ(pnt1.x,pnt1.y,pnt1.z) );

                    // get points on line segment
                    typename pcl::PointCloud<MyPointT>::Ptr on_line_cloud( new pcl::PointCloud<MyPointT>() );
                    {
                        // allocate line
                        pcl::ModelCoefficients::Ptr line_coeffs( new pcl::ModelCoefficients() );
                        line_coeffs->values.resize( 6 );

                        // get starting point ( pnt0 = cloud.back() for vert_id == 0, otherwise pnt0 = cloud[vert_id-1] )
                        MyPointT pnt0 = cloud_hull->at( polygons[poly_id].vertices[vert_id ? vert_id-1 : polygons[poly_id].vertices.size()-1] );
                        std::copy( pnt0.data, pnt0.data+3, line_coeffs->values.begin() ); // put into line primitive

                        // get line direction
                        Eigen::Vector3f direction( pnt1.x - pnt0.x, pnt1.y - pnt0.y, pnt1.z - pnt0.z );
                        float line_length = direction.norm();   // save length for later
                        direction /= line_length;               // normalize direction
                        std::copy( direction.data(), direction.data() + 3, line_coeffs->values.begin()+3 ); // put into line primitive

                        // select line inliers
                        pcl::PointIndices::Ptr line_inliers( new pcl::PointIndices() );
                        pcl::SampleConsensusModelLine<MyPointT> sac_line( cloud_projected );
                        sac_line.selectWithinDistance( coeffs2Vector(line_coeffs), /*inlier distance: */ .005, line_inliers->indices );

                        // project found inlier points to line
                        pcl::ProjectInliers<MyPointT> proj;
                        proj.setInputCloud( cloud_projected );
                        proj.setIndices( line_inliers );
                        proj.setModelCoefficients( line_coeffs );
                        proj.setModelType( pcl::SACMODEL_LINE );
                        proj.filter( *on_line_cloud );

                        // copy points, that are close to hull_line and their projections fall into the line segment
                        for ( size_t point_id = 0; point_id != on_line_cloud->size(); ++point_id )
                        {
                            // check if pnt - pnt0 < pnt1-pnt0  ~  if the projected point is inside the line segment
                            // this is needed due to the concavity of the hull. Because otherwise, this point belongs
                            // to a different hull segment, that it might not be close enough to.
                            float projection = (on_line_cloud->at(point_id).getVector3fMap() - pnt0.getVector3fMap()).dot(direction);
                            if ( (projection <= line_length) && (projection >= 0.f) )
                            {
                                // add original point to output
                                hull_cloud->push_back( cloud_projected->at(line_inliers->indices[point_id]) );
                            }
                        } // for points
                    }
                }  // for vertices

                // add generating polygon to both viewports
                char title[255];
                sprintf( title, "polygon%04lu", poly_id );
                vptr->addPolygon<pcl::PointXYZ>( polygon_cloud_ptr, 1., .2, .2, title, vp0 );
                vptr->addPolygon<pcl::PointXYZ>( polygon_cloud_ptr, 1., .2, .2, title, vp1 );
                //vptr->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_SURFACE, std::string(title) );
            } // for polygons
        } // qhull

        // show qhlull
        vptr->addPointCloud( hull_cloud, "hull", vp1 );
        vptr->spin();

        // get local neighs
        pcl::PointIndices::Ptr indices = smartgeometry::allIndicesOf( hull_cloud );
        std::vector< std::vector<int> > neighbour_indices;
        //getNeighbourhoodIndices( neighbour_indices, hull_cloud, indices );
        smartgeometry::getNeighbourhoodIndices( neighbour_indices
                                                , hull_cloud
                                                , &(indices->indices)
                                                , /* sqr distances: */ NULL
                                                , /*             K: */ 15
                                                , /*    max_radius: */ .1f );

        // get line models using line ransac
        std::vector<pcl::ModelCoefficients::Ptr> models;
        getLocalModels( models
                        , neighbour_indices
                        , hull_cloud
                        );

        // display
        Eigen::Vector3f plane_normal( plane_coeffs->values[0], plane_coeffs->values[1], plane_coeffs->values[2] );
        assert( models.size() == indices->indices.size() );
        for ( size_t i = 0; i < models.size(); i+=3 )
        {
            MyPointT pnt = hull_cloud->at( indices->indices[i] );
            Eigen::Vector3f line_dir( models[i]->values[3], models[i]->values[4], models[i]->values[5] );
            Eigen::Vector3f normal = -1.f*line_dir.cross(plane_normal).normalized();
            char tit[255];
            sprintf( tit, "line_normal%d", static_cast<int>(i) );
            vptr->addArrow( smartgeometry::toPointXYZ(pnt.getVector3fMap() + normal * 0.05f)
                            , pnt
                            , 1., .1, .1
                            , false
                            , tit
                            , vp1 );
        }
        vptr->spin();

        return EXIT_SUCCESS;
    }
} // ns am
#endif

#endif // CUBOIDRANSAC_H
