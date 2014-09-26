#ifndef __GF2_SCHNABELENV_H__
#define __GF2_SCHNABELENV_H__

#include <vector>
#include "pcl/point_cloud.h"
#include "globfit2/my_types.h"

namespace GF2
{
    //class PlanePrimitive;
    typedef pcl::PointNormal MyPoint;
    typedef pcl::PointCloud<MyPoint> MyCloud;

    class SchnabelEnv
    {
        public:
            template <typename PrimitiveT>
            static inline int
            run( std::vector<PrimitiveT>    &planes
                 , pcl::PointCloud<GF2::MyPoint>::Ptr  cloud
                 , float scale = 0.01
                 , int                                 min_support_arg = 300
                 , int show = 1 );
    };

} // ns am

#endif // SCHNABELENV_H
