#ifndef __RAPTER_PCL_UTIL_HPP__
#define __RAPTER_PCL_UTIL_HPP__

#ifdef RAPTER_USE_PCL

#include "Eigen/Eigenvalues" // SelfAdjointSolver
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

// from pcltools
namespace smartgeometry {

    template <typename PointsT, typename Scalar> inline int
    computeCentroid( Eigen::Matrix<Scalar,4,1>       & centroid
                     , PointsT                  const& cloud
                     , std::vector<int>         const* indices_arg )
    {
        centroid.setZero();

        const size_t N = indices_arg ? indices_arg->size() : cloud.size();
        for (size_t pid = 0; pid != N; ++pid )
        {
            const size_t index = indices_arg ? (*indices_arg)[pid] : pid;
            centroid[0] += cloud[index].x;
            centroid[1] += cloud[index].y;
            centroid[2] += cloud[index].z;
        }
        centroid /= static_cast<Scalar>( N );

        return EXIT_SUCCESS;
    } //...computeCentroid()

    // untested with indices
    template <typename PointsT, typename Scalar> inline int
    computeCovarianceMatrix( Eigen::Matrix<Scalar, 3, 3>         & covariance_matrix
                             , PointsT                      const& cloud
                             , std::vector<int>             const* indices_arg
                             , Eigen::Matrix<Scalar, 4, 1>  const* centroid_arg
                             , std::vector<Scalar>          const* weights_arg
                             )
    {
        // Initialize to 0
        covariance_matrix.setZero();

        const size_t N = indices_arg ? indices_arg->size() : cloud.size();

        // init centroid
        Eigen::Matrix<Scalar,4,1> centroid; centroid.setZero();
        if ( centroid_arg )
            centroid = *centroid_arg;
        else
            computeCentroid( centroid, cloud, indices_arg );

        // For each point in the cloud
        for ( size_t pid = 0; pid != N; ++pid )
        {
            const int    index  = indices_arg ? (*indices_arg)[pid  ] : pid;
            const Scalar weight = weights_arg ? (*weights_arg)[pid  ] : 1;

            Eigen::Matrix<Scalar, 4, 1> pt;
            pt[0] = cloud[index].x - centroid[0];
            pt[1] = cloud[index].y - centroid[1];
            pt[2] = cloud[index].z - centroid[2];

            covariance_matrix (1, 1) += weight * pt.y () * pt.y ();
            covariance_matrix (1, 2) += weight * pt.y () * pt.z ();
            covariance_matrix (2, 2) += weight * pt.z () * pt.z ();

            pt *= pt.x ();

            covariance_matrix (0, 0) += weight * pt.x ();
            covariance_matrix (0, 1) += weight * pt.y ();
            covariance_matrix (0, 2) += weight * pt.z ();
        }
        covariance_matrix (1, 0) = covariance_matrix( 0, 1);
        covariance_matrix (2, 0) = covariance_matrix( 0, 2);
        covariance_matrix (2, 1) = covariance_matrix( 1, 2);

        if ( weights_arg )
        {
            Scalar sum_weight = std::accumulate(weights_arg->begin(), weights_arg->end(), 0.);
            if ( sum_weight > FLT_EPSILON ) covariance_matrix /= sum_weight;
        }
        else
        {
            if ( N == 0 ) std::cerr << "[" << __func__ << "]: " << "dividing by N: " << N << std::endl;
            else covariance_matrix /= static_cast<Scalar>( N );
        }

        return EXIT_SUCCESS;
    } //...computeCovarianceMatrix()

    namespace geometry
    {
        // Point2Primitive distance
        template <typename Scalar, int rows> Scalar
        pointPrimitiveDistance (Eigen::Matrix<Scalar,3,1> const& pnt,Eigen::Matrix<Scalar,rows,1> const& primitive);
        template<> inline float
        pointPrimitiveDistance<float,6> (Eigen::Matrix<float,3,1> const& pnt, Eigen::Matrix<float,6,1> const& line )
        {
            return (line.template head<3>() - pnt).cross( line.template segment<3>(3) ).norm();
        }
        template<> inline float
        pointPrimitiveDistance<float,4> (Eigen::Matrix<float,3,1> const& pnt, Eigen::Matrix<float,4,1> const& plane )
        {
            return plane.template head<3>().dot( pnt ) + plane(3);
        }

        // Primitive from point and normal
        template <typename Scalar, int rows> Eigen::Matrix<Scalar,rows,1>
        fromPointAndNormal( Eigen::Matrix<Scalar,3,1> const& pnt,Eigen::Matrix<Scalar,3,1> const& normal );
        template <> inline Eigen::Matrix<float,6,1>
        fromPointAndNormal<float,6>( Eigen::Matrix<float,3,1> const& pnt,Eigen::Matrix<float,3,1> const& normal )
        {
            return (Eigen::Matrix<float,6,1>() << pnt, normal).finished(); // TODO: this is bullshit, it's not the normal, but the direction...
        }
        template <> inline Eigen::Matrix<float,4,1>
        fromPointAndNormal<float,4>( Eigen::Matrix<float,3,1> const& pnt,Eigen::Matrix<float,3,1> const& normal )
        {
            //model_coefficients[3] = -1 * (model_coefficients.template head<4>().dot (p0.matrix ()));
            Eigen::Matrix<float,4,1> primitive;
            primitive.template segment<3>(0) = normal;
            primitive                    (3) = static_cast<float>(-1) * primitive.template head<3>().dot( pnt.template head<3>() ); // distance

            return primitive;
        }

        /**
         * @brief fitLine               [Re]Fits 3D line to a [part of a] pointcloud.
         * @param line                  Output line, and possibly input line to refit, if \param start_from_input_line is true.
         * @param cloud                 Points to fit to. Must have methods operator[] and getVector3fMap()->Eigen::Vector3f. If \param p_indices!=NULL, must have at least max(*p_indices) points.
         * @param scale                 Distance where point get's zero weight
         * @param p_indices             Indices to use from cloud. Can be NULL, in which case the whole cloud is used.
         * @param refit                 How many refit iterations. 0 means once, obviously (TODO to fix...).
         * @param start_from_input_line Assume, that \param line contains a meaningful input, and calculate weights on the 0th iteration already.
         */
        template <class PointsT, typename Scalar = float, int rows = 6, class _DerivedT = Eigen::Matrix<Scalar,4,1> > inline int
        fitLinearPrimitive( _DerivedT    & primitive // Eigen::Matrix<Scalar,rows,1>
                            , PointsT                  const& cloud
                            , Scalar                          scale
                            , std::vector<int>              * p_indices             = NULL
                            , int                             refit                 = 0
                            , bool                            start_from_input      = false
                            , Scalar                       (*pointPrimitiveDistanceFunc)(Eigen::Matrix<Scalar,3,1> const& pnt, Eigen::Matrix<Scalar,rows,1> const& primitive) = &(pointPrimitiveDistance<Scalar,rows>)
                            , Eigen::Matrix<Scalar,rows,1> (*    fromPointAndNormalFunc)(Eigen::Matrix<Scalar,3,1> const& pnt, Eigen::Matrix<Scalar,3   ,1> const& normal   ) = &(fromPointAndNormal<Scalar,rows>)
                            , bool                            debug                 = false )
        {
            eigen_assert( (rows == 4) || (rows == 6) );

            // number of points to take into account
            const size_t N = p_indices ? p_indices->size() : cloud.size();

            // skip, if not enought points found to fit to
            if ( N < 2 ) { std::cerr << "[" << __func__ << "]: " << "can't fit line to less then 2 points..." << std::endl; return EXIT_FAILURE; }

            int iteration = 0; // track refit iterations
            do
            {
                // LeastSquares weights
                std::vector<Scalar> weights( N, 1.f );

                // calculate weights, if value in "line" already meaningful
                if ( start_from_input || (iteration > 0) )
                {
                    // calculate distance from all points
                    for ( size_t point_id = 0; point_id != N; ++point_id )
                    {
                        // formula borrowed from PCL: (line_pt - point).cross3(line_dir).squaredNorm();
                        Eigen::Matrix<Scalar,3,1> pnt = cloud[ p_indices ? (*p_indices)[point_id] : point_id ].getVector3fMap();
                        weights[point_id] = pointPrimitiveDistanceFunc( pnt, primitive );
                    }

                    // the farther away, the smaller weight -->
                    // w_i = f( dist_i / scale ), dist_i < scale; f(x) = (x^2-1)^2
                    for ( size_t wi = 0; wi != weights.size(); ++wi )
                    {
                        if ( weights[wi] < scale )
                        {
                            weights[wi] /= scale;                                               // x = dist_i / scale
                            weights[wi] = (weights[wi] * weights[wi] - static_cast<Scalar>(1)); // x^2-1
                            weights[wi] *= weights[wi];                                         // (x^2-1)^2
                        }
                        else
                            weights[wi] = static_cast<Scalar>(0);                               // outside scale, truncated to 0
                    }
                }

                // compute centroid of cloud or selected points
                Eigen::Matrix<Scalar,4,1> centroid;
                smartgeometry::computeCentroid( centroid, cloud, p_indices );

                // compute neighbourhood covariance matrix
                Eigen::Matrix<Scalar,3,3> cov;
                smartgeometry::computeCovarianceMatrix( cov, cloud, p_indices, &centroid, &weights ); // weights might be all 1-s

                // solve for neighbourhood biggest eigen value
                Eigen::SelfAdjointEigenSolver< Eigen::Matrix<Scalar, 3, 3> > es;
                es.compute( cov );

                if ( rows == 6 ) // line -> dir ==
                {
                    // get eigen vector for biggest eigen value
                    const int max_eig_val_id = std::distance( es.eigenvalues().data(), std::max_element( es.eigenvalues().data(), es.eigenvalues().data()+3 ) );

                    // output line
                    primitive = fromPointAndNormalFunc( centroid.template head<3>(),
                                                        es.eigenvectors().col(max_eig_val_id).normalized() );
                }
                else if ( rows == 4 ) // plane
                {
                    // get eigen vector for biggest eigen value
                    const int min_eig_val_id = std::distance( es.eigenvalues().data(), std::min_element( es.eigenvalues().data(), es.eigenvalues().data()+3 ) );

                    // output line
                    primitive = fromPointAndNormalFunc( centroid.template head<3>(),
                                                        es.eigenvectors().col(min_eig_val_id).normalized() );
                }
                else
                    std::cerr << "[" << __func__ << "]: " << "lines(rows==6) or planes(rows==4), not rows == " << rows << std::endl;

            }
            while ( iteration++ < refit );

            return EXIT_SUCCESS;
        } // ... fitline
    } // ... ns geometry
} // ...ns smartgeometry

#endif // RAPTER_USE_PCL

#endif // __RAPTER_PCL_UTIL_HPP__
