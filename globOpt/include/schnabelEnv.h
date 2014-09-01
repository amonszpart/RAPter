#ifndef __GF2_SCHNABELENV_H__
#define __GF2_SCHNABELENV_H__

#include <vector>
#include "pcl/point_cloud.h"
#include "my_types.h"

namespace am
{
    class PlanePrimitive;

    class SchnabelEnv
    {
        public:
            static int
            run( std::vector<::am::PlanePrimitive>    &planes
                 , pcl::PointCloud<MyPoint>::ConstPtr  cloud
                 , int                                 min_support_arg = 300 );
    };

} // ns am

#endif // SCHNABELENV_H
