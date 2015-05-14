#ifndef __RAPTER_PCL_UTIL_HPP__
#define __RAPTER_PCL_UTIL_HPP__

#ifdef RAPTER_USE_PCL

#include "pcl/point_types.h"
#include "pcl/point_cloud.h"
#include "pcl/search/kdtree.h"

namespace rapter {
    namespace pclutil {

        typedef          pcl::PointXYZ                          PclSearchPointT;
        typedef          pcl::search::KdTree<PclSearchPointT>   PclSearchTreeT;
        typedef typename PclSearchTreeT::Ptr                    PclSearchTreePtrT;

        template <class PointAllocatorFunctorT, class PointContainerT, class CloudPtrT> inline int
        cloudToVector( PointContainerT &container, CloudPtrT const& cloud )
        {
            const size_t count = cloud->size();
            container.reserve( count );
            for ( size_t pid = 0; pid != count; ++pid )
            {
                container.push_back( PointAllocatorFunctorT::eval(cloud->at(pid).getVector3fMap()) );
            }

            return EXIT_SUCCESS;
        } //...cloudToVector

        template <typename Derived>
        inline ::pcl::PointXYZ
        asPointXYZ( Derived const& vector3 )
        {
            return ::pcl::PointXYZ( vector3.x(), vector3.y(), vector3.z() );
        }

        template <class _PointContainerT>
        inline PclSearchTreePtrT buildANN( _PointContainerT const& points )
        {
            pcl::PointCloud<PclSearchPointT>::Ptr ann_cloud( new pcl::PointCloud<PclSearchPointT>() );
            {
                ann_cloud->reserve( points.size() );
                for ( size_t pid = 0; pid != points.size(); ++pid )
                {
                    PclSearchPointT pnt;
                    pnt.x = points[pid].template pos()(0);
                    pnt.y = points[pid].template pos()(1);
                    pnt.z = points[pid].template pos()(2);
                    ann_cloud->push_back( pnt );
                }
            }

            PclSearchTreePtrT tree( new PclSearchTreeT() );
            tree->setInputCloud( ann_cloud );

            return tree;
        } //...buildANN

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

    } //...ns pclutil
} //...ns rapter

#endif // RAPTER_USE_PCL

#endif // __RAPTER_PCL_UTIL_HPP__
