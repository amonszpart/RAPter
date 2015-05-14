#ifndef __RAPTER_PCLUTIL_H__
#define __RAPTER_PCLUTIL_H__

#ifdef RAPTER_USE_PCL

#include "pcl/point_types.h"
#include "pcl/point_cloud.h"

namespace rapter {
    typedef          pcl::PointNormal                PclPointT;
    typedef          pcl::PointCloud<PclPointT>      PclCloudT;
    typedef typename pcl::PointCloud<PclPointT>::Ptr PclCloudPtrT;
}

#endif // RAPTER_USE_PCL

#endif // __RAPTER_PCLUTIL_H__
