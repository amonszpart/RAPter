#ifndef __GF2_GTCREATOR_H__
#define __GF2_GTCREATOR_H__

#include <string>
#include "pcl/point_cloud.h"
#include "globfit2/my_types.h"
#include <opencv2/core.hpp>

namespace GF2
{

    class GTCreator
    {
    public:
        typedef pcl::PointXYZRGB PointT;

        int
        image_2_2DCloud( pcl::PointCloud<PointT>::Ptr &cloud
                         , cv::Mat                const& img
                         , int                         N_samples
                         , const float                 Z          );

        //static int run( pcl::PointCloud<PointT>::Ptr &cloud, std::string gt_name, int gt_nPoints, float gt_noise );
        static int sampleImage( pcl::PointCloud<PointT>::Ptr &cloud
                                , std::string image_path
                                , int gt_nPoints
                                , float gt_noise
                                , float gt_scene_size            = 1.f
                , float filter_size              = 0.075f
                , Eigen::Vector3f *sensor_origin = NULL
                );
    };

} // ns GF2

#endif // __GF2_GTCREATOR_H__
