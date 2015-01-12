#ifndef MY_TYPES_H
#define MY_TYPES_H

#include <Eigen/Dense>
#include <set>

#if GF2_USE_PCL
#   include "pcl/point_types.h"
#   include "pcl/point_cloud.h"
#endif

namespace GF2
{
    template<typename Scalar, int Dim> inline Scalar
    angleInRad( Eigen::Matrix<Scalar,Dim,1> const& v1, Eigen::Matrix<Scalar,Dim,1> const& v2 )
    {
        // better, than acos
        Scalar angle = atan2( v1.cross(v2).norm(), v1.dot(v2) );

        // fix nans
        if ( angle != angle )   angle = Scalar(0);

        return angle;
    }
#if 0
    template<typename Scalar, int Dim> inline Scalar
    angleInRadSigned( Eigen::Matrix<Scalar,Dim,1> const& v1, Eigen::Matrix<Scalar,Dim,1> const& v2 )
    {
        // better, than acos
        Scalar angle = atan2( v1.cross(v2).norm(), v1.dot(v2) );

        // fix nans
        if ( angle != angle )   angle = Scalar(0);

        Scalar sign = Eigen::Matrix<Scalar,3,1>::UnitZ()

        return angle;
    }
#endif

#if GF2_USE_PCL

    typedef          pcl::PointNormal                PclPointT;
    typedef          pcl::PointCloud<PclPointT>      PclCloudT;
    typedef typename pcl::PointCloud<PclPointT>::Ptr PclCloudPtrT;

    template <int Dim>
    struct PCLPointAllocator
    {
            template <class PCLPointT, class FromPointT>
            PCLPointT
            static inline create( FromPointT const& pnt );
    };

    template <>
    template <class PCLPointT, class FromtPointT>
    PCLPointT
    PCLPointAllocator<3>::create( FromtPointT const& pnt )
    {
        PCLPointT pcl_pnt;
        pcl_pnt.x = pnt(0);
        pcl_pnt.y = pnt(1);
        pcl_pnt.z = pnt(2);

        return pcl_pnt;
    }

    template <>
    template <class PCLPointT, class FromtPointT>
    PCLPointT
    PCLPointAllocator<6>::create( FromtPointT const& pnt ) // x, y, z, nx, ny, nz
    {
        PCLPointT pcl_pnt;
        pcl_pnt.x = pnt(0);
        pcl_pnt.y = pnt(1);
        pcl_pnt.z = pnt(2);

        return pcl_pnt;
    }

#endif
}

#endif // MY_TYPES_H
