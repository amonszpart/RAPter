#include "cuboidRansac.h"

#include "my_types.h"
#include "pcl/point_cloud.h"
#include "pcl/point_types.h"
#include "pcl/PointIndices.h"
#include "pcl/io/pcd_io.h"
#include "pcl/filters/project_inliers.h"
#include "pcltools/primitives.h" // smartgeom
#include "pcl/segmentation/sac_segmentation.h"
#include "pcl/PointIndices.h"
#include "pcl/search/kdtree.h"

#include "pcl/surface/concave_hull.h"

namespace am
{

    int
    CuboidRansac::getNeighbourhoodIndices( std::vector<pcl::PointIndices::Ptr>  &neighbour_indices
                                           , pcl::PointCloud<MyPoint>::ConstPtr  cloud
                                           , pcl::PointIndices::ConstPtr         indices_arg
                                           , std::vector<std::vector<float> >   *distances
                                           , int                                 K                  )
    {
        // prepare output
        pcl::PointIndices::ConstPtr indices = indices_arg ? indices_arg : smartgeometry::allIndicesOf( cloud );
        neighbour_indices.resize( indices->indices.size() );
        if ( distances )
            distances->resize( indices->indices.size() );

        pcl::search::KdTree<MyPoint>::Ptr tree( new pcl::search::KdTree<MyPoint> );
        pcl::IndicesPtr indices_ptr( new std::vector<int>() );
        *indices_ptr = indices->indices;
        tree->setInputCloud( cloud, indices_ptr );

        std::vector<int>   neighs   ( K );
        std::vector<float> sqr_dists( K );

        for ( size_t i = 0; i != indices->indices.size(); ++i )
        {
            MyPoint searchPoint = cloud->at( indices->indices[i] );
            if ( tree->nearestKSearch(searchPoint, K, neighs, sqr_dists) > 0 )
            {
                neighbour_indices[i] = pcl::PointIndices::Ptr( new pcl::PointIndices() );
                neighbour_indices[i]->indices = neighs;

                if ( distances )
                    distances->at(i) = sqr_dists;

//                for ( size_t i = 0; i < neighs.size (); ++i )
//                    std::cout << "    "  <<   cloud->points[ neighs[i] ].x
//                              << " " << cloud->points[ neighs[i] ].y
//                              << " " << cloud->points[ neighs[i] ].z
//                              << " (squared distance: " << sqr_dists[i] << ")" << std::endl;
            }
            else
                std::cout << __func__ << ": no neighs found" << std::endl;
        }

        return EXIT_SUCCESS;
    }

    int
    CuboidRansac::getLocalModels( std::vector<pcl::ModelCoefficients::Ptr>    & models
                                 , std::vector<std::vector<int> >       const & neighbour_indices
                                 , pcl::PointCloud<MyPoint>::ConstPtr           cloud
                                 , int                                          sac_model
                                 , double                                       distanceThreshold )
    {
        models.reserve( neighbour_indices.size() );

        pcl::SACSegmentation<MyPoint> seg;
        seg.setModelType        ( sac_model          );
        seg.setDistanceThreshold( distanceThreshold  );
        seg.setMethodType       ( pcl::SAC_RANSAC    );
        seg.setInputCloud       ( cloud              );

        pcl::PointIndices inliers;
        pcl::PointIndices::Ptr neighbours( new pcl::PointIndices() );
        for ( size_t i = 0; i != neighbour_indices.size(); ++i )
        {
            models.emplace_back( pcl::ModelCoefficients::Ptr(new pcl::ModelCoefficients()) );
            neighbours->indices = neighbour_indices[i];
            seg.setIndices( neighbours );
            seg.segment   ( inliers, *models[i] );
        }

        return EXIT_SUCCESS;
    }

#if 0
    template <typename MyPointT>
    int Segmentation::extract_plane( pcl::ModelCoefficients::Ptr                        &out_coefficients,
                                     boost::shared_ptr<const pcl::PointCloud<MyPointT>> in_cloud_ptr,
                                     pcl::PointIndices::ConstPtr                        in_indices_ptr,
                                     double                                             distanceThreshold,
                                     int                                                sac_model,
                                     Eigen::Vector3f                                    *p_normal            )
    {
        // prepare output
        if ( !out_coefficients )
            out_coefficients = boost::make_shared< pcl::ModelCoefficients >();

        // Normals
        pcl::PointCloud<pcl::Normal>::Ptr roi_cloud_normals_ptr (new pcl::PointCloud<pcl::Normal>);
        {
            pcl::NormalEstimation<MyPointT, pcl::Normal> normal_estimation;
            typename pcl::search::KdTree<MyPointT>::Ptr  tree( new pcl::search::KdTree<MyPointT> );

            normal_estimation.setSearchMethod( tree           );
            normal_estimation.setInputCloud  ( in_cloud_ptr   );
            normal_estimation.setKSearch     ( 20            );
            if ( in_indices_ptr )
                normal_estimation.setIndices( in_indices_ptr );

            normal_estimation.compute        ( *roi_cloud_normals_ptr );
        }

        // resize normals cloud to match cloud_ptr, copy estimated ones to right positions, leave others as 0
        pcl::PointCloud<pcl::Normal>::Ptr cloud_normals_ptr ( new pcl::PointCloud<pcl::Normal> );
        if ( in_indices_ptr )
        {
            _expandNormals( /*   out_cloud: */ cloud_normals_ptr,
                            /*   in_indice: */ in_indices_ptr,
                            /*  in_normals: */ roi_cloud_normals_ptr,
                            /* output_size: */ in_cloud_ptr->size() );
        }

        // Segmentation
        pcl::PointIndices::Ptr inliers_indices( new pcl::PointIndices );
        pcl::SACSegmentationFromNormals<MyPointT, pcl::Normal> seg;
        // in
        seg.setInputCloud           ( in_cloud_ptr                   );
        seg.setInputNormals         ( in_indices_ptr ? cloud_normals_ptr : roi_cloud_normals_ptr );
        if ( in_indices_ptr ) seg.setIndices( in_indices_ptr );
        // props
        seg.setOptimizeCoefficients ( true                       );
        seg.setMethodType           ( pcl::SAC_RANSAC            );
        seg.setDistanceThreshold    ( distanceThreshold          ); // .1, .05
        seg.setMaxIterations        ( 500                        );

        switch ( sac_model )
        {
            case pcl::SACMODEL_PERPENDICULAR_PLANE:
                seg.setModelType ( pcl::SACMODEL_PARALLEL_PLANE );
                seg.setAxis      ( *p_normal                    );
                seg.setEpsAngle  ( pcl::deg2rad(20.f)           );
                break;
            case pcl::SACMODEL_NORMAL_PLANE:
                seg.setModelType ( pcl::SACMODEL_NORMAL_PLANE );
                seg.setNormalDistanceWeight ( 0.1                        );
                break;
            case pcl::SACMODEL_NORMAL_PARALLEL_PLANE:
                seg.setModelType ( pcl::SACMODEL_NORMAL_PARALLEL_PLANE );
                seg.setAxis      ( *p_normal                    );
                seg.setEpsAngle  ( pcl::deg2rad(5.f)           );
                break;
            case pcl::SACMODEL_PARALLEL_PLANE:
                seg.setModelType ( pcl::SACMODEL_PARALLEL_PLANE );
                seg.setAxis      ( *p_normal                    );
                seg.setEpsAngle  ( pcl::deg2rad(5.f)           );
                break;
            default:
                std::cerr << "sac_model_plane unknown" << std::endl;
                break;
        }

        // out
        seg.segment ( *inliers_indices, *out_coefficients );
        return (out_coefficients->values.size() > 0) ? EXIT_SUCCESS : EXIT_FAILURE;
    }
#endif

}
