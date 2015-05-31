#ifndef __RAPTER_SCHNABELENV_H__
#define __RAPTER_SCHNABELENV_H__

#include <vector>
#include "pcl/point_cloud.h"
#include "rapter/typedefs.h"

namespace rapter
{
    //class PlanePrimitive;

    class SchnabelEnv
    {
        public:
            template <class PclCloudT, typename PrimitiveT, /*class PidGidT, */class PointContainerT >
            static inline int
            run( std::vector<PrimitiveT>    & planes
                 //, PidGidT                & pidGid
                 , PointContainerT          & outPoints
                 , PointContainerT     const& points
                 , typename PclCloudT::Ptr  & cloud
                 , float scale = 0.01
                 , int                                 min_support_arg = 300
                 , int show = 1
                 , bool extrude2D = false
                 , int pointMultiplier = 50
                 );
    };

} // ns am

#endif // __RAPTER_SCHNABELENV_H__
