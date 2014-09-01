#ifndef __GF2_PLANEPRIMITIVE_H__
#define __GF2_PLANEPRIMITIVE_H__

#include <Eigen/Dense>
#include "primitives/primitive.h"
#include "pcltools/util.hpp"

#ifdef GF2_USE_PCL
#   include "pcl/ModelCoefficients.h"
#   include "pcl/sample_consensus/sac_model_plane.h"
#   include "pcl/common/common.h"
#   include "pcl/visualization/pcl_visualizer.h"
#endif

namespace GF2
{
    //! \brief   Class to wrap a plane with. Implements #pos() and #dir() functions, and a constructor taking a position and a direction.
    //!
    //!          Stores 3D normal at the first three coeffs, and distance from origin at the fourth coordinate.
    class PlanePrimitive : public ::GF2::Primitive<4>
    {
        public:
            //! \brief Inherited constructor from Primitive.
            using ::GF2::Primitive<4>::Primitive;

            // ____________________CONSTRUCT____________________

            //! \brief           Creates PlanePrimitive from point on plane and direction.
            //! \param[in] p0    Point on plane.
            //! \param[in] dir   Plane normal.
            PlanePrimitive( Eigen::Matrix<Scalar,3,1> pnt, Eigen::Matrix<Scalar,3,1> normal );

            // ____________________VIRTUALS____________________
            //! \brief  Compulsory virtual overload of position getter. The position of the plane is calculated on the fly from the formula N . x0 + d = 0.
            //! \return The position of the plane as a 3D Eigen::Vector.
            virtual Eigen::Matrix<Scalar,3,1> pos() const { return _coeffs.head<3>() * (-1.f*_coeffs(3)); }
            //! \brief  Compulsory virtual overload of orientation getter. The orientation of the plane is the normal stored at the first three coordinates of #_coeffs.
            //! \return The normal of the plane as a 3D Eigen::Vector.
            virtual Eigen::Matrix<Scalar,3,1> dir() const { return _coeffs.head<3>()                    ; }

            //! \brief  Returns the normal, that is stored at the first three coordinates of the internal storage.
            //! \return The plane normal as a 3D Eigen::Vector Map.
            inline Eigen::Map<Eigen::Matrix<Scalar,3,1> > const normal() { return _coeffs.data(); }

            //! \brief              Returns point to plane distance.
            //! \param[in] point    Point to calculate distance from.
            //! \return             Distance from point to plane.
            inline Scalar
            getDistance( Eigen::Matrix<Scalar,3,1> const& point )
            {
                return this->dir().dot(point) + this->_coeffs(3);
            }

            //! \brief                          Calculates the length of the plane based on the points in \p cloud, masked by \p indices and the distance from point to plane \p threshold.
            //!
            //!                                 The method calculates the inliers and selects the most far away point from #pos() in both directions.
            //! \tparam     _PointT             Type of endpoints returned in \p minMax. Concept: _PointContainerT::element_type::PointType == pcl::PointXYZRGB.
            //! \tparam     _PointContainerPtrT Type to store the points to select inliers from. Concept: pcl::PointCloud< _PointT >::Ptr.
            //! \param[out] minMax              Output container holding the four endpoints.
            //! \param[in]  cloud               Point container to look for inlier points in.
            //! \param[in]  threshold           A point in \p cloud is an inlier, if it is closer than this threshold. Usually the "scale".
            //! \param[in]  indices             Optional input to specify subset of points by indices.
            //! \return                         EXIT_SUCCESS
            //! \todo                           Detach from PCL.
            template <typename _PointT, typename _PointContainerPtrT>
            int
            getExtent( std::vector<_PointT>        & minMax
                       , _PointContainerPtrT         cloud
                       , const double                threshold = 0.01
                       , std::vector<int>     const* indices   = NULL ) const
            {
#ifdef GF2_USE_PCL
                pcl::PointIndices::Ptr                inliers( new pcl::PointIndices() );
                typename pcl::PointCloud<_PointT>::Ptr on_plane_cloud( new pcl::PointCloud<_PointT>() );
                {
                    pcl::SampleConsensusModelPlane<_PointT> sacplane( cloud );
                    if ( indices )
                        sacplane.setIndices( *indices );
                    sacplane.selectWithinDistance(          coeffs(), threshold, inliers->indices      );
                    sacplane.projectPoints       ( inliers->indices ,  coeffs(), *on_plane_cloud, false );
                }

                if ( !on_plane_cloud->size() )
                {
                    return EXIT_FAILURE;
                }

                Eigen::Matrix<float,4,4> frame; // 3 major vectors as columns, and the fourth is the centroid
                {
                    smartgeometry::PCA( frame, on_plane_cloud, NULL );

                    // get the unit axis that is most perpendicular to the 3rd dimension of the frame
                    std::pair<Eigen::Vector3f,float> dim3 = { Eigen::Vector3f::Zero(), FLT_MAX };
                    {
                        float tmp;
                        if ( (tmp=std::abs(frame.col(2).head<3>().dot( Eigen::Vector3f::UnitX() ))) < dim3.second ) { dim3.first = Eigen::Vector3f::UnitX(); dim3.second = tmp; }
                        if ( (tmp=std::abs(frame.col(2).head<3>().dot( Eigen::Vector3f::UnitY() ))) < dim3.second ) { dim3.first = Eigen::Vector3f::UnitY(); dim3.second = tmp; }
                        if ( (tmp=std::abs(frame.col(2).head<3>().dot( Eigen::Vector3f::UnitZ() ))) < dim3.second ) { dim3.first = Eigen::Vector3f::UnitZ(); dim3.second = tmp; }
                    }
                    frame.col(0).head<3>() = frame.col(2).head<3>().cross( dim3.first ).normalized();
                    frame.col(1).head<3>() = frame.col(2).head<3>().cross( frame.col(0).head<3>() ).normalized();
                }

                typename pcl::PointCloud<_PointT>::Ptr local_cloud = NULL;
                smartgeometry::cloud2Local<_PointT>( local_cloud, frame, on_plane_cloud, pcl::PointIndices::ConstPtr() );

                _PointT min_pt, max_pt;
                pcl::getMinMax3D( *local_cloud, min_pt, max_pt );

                minMax.resize( 4 );
                minMax[0] = minMax[1] = min_pt;
                minMax[1].y = max_pt.y;
                minMax[2] = minMax[3] = max_pt;
                minMax[3].y = min_pt.y;

                for ( int d = 0; d != 4; ++d )
                {
                    // to world
                    Eigen::Matrix<float,4,1> pnt = frame * minMax[d].getVector4fMap();
                    // copy
                    minMax[d].x = pnt(0); minMax[d].y = pnt(1); minMax[d].z = pnt(2);
                }

                return EXIT_SUCCESS;
#           else
                std::cerr << "[" << __func__ << "]: " << "Needs PCL to work!!!!"
                return exit_FAILURE;
#           endif //...GF2_USE_PCL
            } //...getExtent()

            // ____________________DRAWING____________________
#ifdef GF2_USE_PCL
            inline pcl::ModelCoefficients::Ptr
            modelCoefficients() const;

            template <class PointsT> static inline int
            draw( PointsT const& corners
                  , pcl::visualization::PCLVisualizer::Ptr v
                  , std::string plane_name
                  , double r, double g, double b
                  , int viewport_id = 0
                  )
            {
                if ( corners.size() != 4 )
                {
                    std::cerr << "[" << __func__ << "]: " << "...cannot draw with not 4 corners" << std::endl;
                    return EXIT_FAILURE;
                }
                Eigen::Vector3f centroid( Eigen::Vector3f::Zero() );
                for ( size_t corner_id = 0; corner_id != corners.size(); ++corner_id )
                {
                    centroid += corners[corner_id].getVector3fMap();
                }
                centroid /= static_cast<float>(corners.size());

                // polygon containing the for corners spanning the plane's face
                pcl::PointCloud<pcl::PointXYZ>::Ptr plane_polygon_cloud_ptr( new pcl::PointCloud<pcl::PointXYZ> );

                //char title[255];
                for ( int corner_id = 0;
                      corner_id < corners.size();
                      ++corner_id )
                {
                    plane_polygon_cloud_ptr->push_back( corners[corner_id] );
                }

                // draw plane polygon
                v->addPolygon<pcl::PointXYZ>( plane_polygon_cloud_ptr, r,g,b, plane_name, viewport_id );
                v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, .5, plane_name);
                v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_SURFACE, plane_name );

                // show normal
//                if ( drawNormal )
//                {
//                    v->addArrow( am::util::pcl::asPointXYZ( centroid + plane.Normal() * .2f ), // end point
//                                 am::util::pcl::asPointXYZ( centroid                        ), // start point
//                                 r*3., g*3., b*3.,
//                                 false,                                                                   // show length
//                                 std::string(title) + "_normal"                                                  // cloud id
//                                 , viewport_id );
//                }

//                // draw coorindate system
//                if ( drawFrame )
//                {
//                    for ( int axis_id = 0; axis_id < 3; ++axis_id )
//                    {
//                        char arrow_name[1024];
//                        sprintf( arrow_name, "_arrow%d", axis_id );

//                        v->addArrow( am::util::pcl::asPointXYZ( plane.Frame().block<3,1>(0,3) + plane.Frame().block<3,1>(0,axis_id).normalized() *.2f ), // end point
//                                     am::util::pcl::asPointXYZ( plane.Frame().block<3,1>(0,3)                             ), // start point
//                                     (2-axis_id  )%3 *.45 + .1,                                               // r
//                                     (2-axis_id+1)%3 *.45 + .1,                                               // g
//                                     (2-axis_id+2)%3 *.45 + .1,                                               // b
//                                     false,                                                                   // show length
//                                     std::string(title) + arrow_name                                          // cloud id
//                                     , viewport_id );
//                        //v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_, 3., std::string(title) + arrow_name );
//                    }
//                }
                return EXIT_SUCCESS;
            }

            template <class PointsPtrT, class PointT = typename PointsPtrT::element_type::PointType, class Scalar = float> static inline int
            draw( PlanePrimitive const& plane
                  , PointsPtrT cloud
                  , Scalar radius
                  , std::vector<int> const* indices
                  , pcl::visualization::PCLVisualizer::Ptr v
                  , std::string plane_name
                  , double r, double g, double b
                  , int viewport_id = 0
                )
            {
                std::vector<PointT> minMax;
                int err = plane.getExtent( minMax
                                          , cloud
                                          , radius
                                          , indices
                                          , /* no_pca: */ true );
                if ( EXIT_SUCCESS != err )
                {
                    v->addPlane( *plane.modelCoefficients(), plane_name, 0 );
                    return err;
                }
                else
                {
                    std::vector<pcl::PointXYZ> ps;
                    for ( int i = 0; i != minMax.size(); ++i )
                    {
                        pcl::PointXYZ pnt;
                        pnt.x = minMax[i].x; pnt.x = minMax[i].y; pnt.y = minMax[i].z;
                        ps.push_back( pnt );
                    }
                    err += draw( ps, v, plane_name, r, g, b, viewport_id );
                }

                return err;
            }

#endif

            // ____________________DEPRECATED____________________
            //! \deprecated Decides, if two planes are different up to a position and and an angle threshold.
            static bool
            different( PlanePrimitive const& me, PlanePrimitive const& other, Scalar pos_diff, Scalar ang_diff )
            {
                if ( fabs(me()(3)-other()(3)) > pos_diff ) return true;

                Scalar angle = acos( me.dir().dot(other.dir()) ); if ( angle != angle ) angle = static_cast<Scalar>(0);
                return ( angle > ang_diff );
            }
#if 0
            static inline Scalar
            point4Distance( Eigen::Matrix<Scalar,4,1> const& point, Eigen::Matrix<Scalar,-1,1> const& plane )
            {
                return plane.template head<3>().dot(point.template head<3>()) + plane(3);
            }
            static inline Scalar
            point3Distance( Eigen::Matrix<Scalar,3,1> const& point, Eigen::Matrix<Scalar,-1,1> const& plane )
            {
                return plane.template head<3>().dot(point.template head<3>()) + plane(3);
            }
            inline Scalar point3Distance(            Eigen::Matrix<Scalar,3,1>   const  point ) const { return point3Distance( point, _coeffs); }
            inline Scalar point4Distance(            Eigen::Matrix<Scalar,4,1>   const  point ) const { return point4Distance( point, _coeffs); }

            // this does not work
            inline Eigen::Matrix<Scalar,6,1>
            toLineCoeffsAtZ( Scalar z );

            inline bool
            intersectWithPlane( Eigen::Matrix<Scalar,-1,1>      & line
                              , Eigen::Matrix<Scalar,4,1>  const& plane
                              , double                            angular_tolerance ) const;

            inline static bool
            intersectWithPlane( Eigen::Matrix<Scalar,-1,1>       & line
                                , Eigen::Matrix<Scalar,4,1> const& plane_a
                                , Eigen::Matrix<Scalar,4,1> const& plane_b
                                , double                           angular_tolerance );
#endif

    }; //...class PlanePrimitive
} //...ns GF2



// ________________________________________________________HPP_________________________________________________________

#ifndef __GF2_INC_PLANEPRIMITIVE_HPP__
#define __GF2_INC_PLANEPRIMITIVE_HPP__
namespace GF2
{
    inline
    PlanePrimitive::PlanePrimitive( Eigen::Matrix<Scalar, 3, 1> pnt, Eigen::Matrix<Scalar, 3, 1> normal )
    {
        _coeffs.template segment<3>(0) = normal.normalized();
        _coeffs                    (3) = static_cast<Scalar>(-1) * _coeffs.template head<3>().dot( pnt.template head<3>() ); // distance
    }

#   ifdef GF2_USE_PCL
    inline pcl::ModelCoefficients::Ptr
    PlanePrimitive::modelCoefficients() const
    {
        pcl::ModelCoefficients::Ptr model_coeffs( new pcl::ModelCoefficients() );
        model_coeffs->values.resize(Dim);
        std::copy( _coeffs.data(), _coeffs.data()+Dim, model_coeffs->values.begin() );

        return model_coeffs;
    }
#   endif // GF2_USE_PCL

#if 0
    // this does not work
    inline Eigen::Matrix<PlanePrimitive::Scalar,6,1>
    PlanePrimitive::toLineCoeffsAtZ( Scalar z )
    {
        Eigen::Matrix<Scalar,6,1> out;
        out.segment<3>( 3 ) = normal();

        Eigen::Matrix<Scalar,4,1> p0_3d; p0_3d << static_cast<Scalar>(1.f),static_cast<Scalar>(1.f),z,static_cast<Scalar>(1.f);
        Scalar dist = p0_3d.dot( _coeffs ); // project point(1,1,z) onto the plane, [1,1,z,1] . [nx,ny,nz,d] = sum([nx,ny,nz,d])
        out.head<3>() = p0_3d.head<3>() - dist * normal();
        out.segment<3>(3).normalize();
        std::cout << "[" << __func__ << "]: " << " plane: " << _coeffs.transpose() << std::endl;
        std::cout << "[" << __func__ << "]: " << " p0 is now " << out.head<3>().transpose() << std::endl;
        std::cout << "[" << __func__ << "]: " << " dir is now " << out.segment<3>(3).transpose() << std::endl;

        return out;
    }
    inline bool
    PlanePrimitive::intersectWithPlane( Eigen::Matrix<Scalar,-1,1>       & line
                                        , Eigen::Matrix<Scalar,4,1> const& plane
                                        , double                           angular_tolerance ) const
    {
        return PlanePrimitive::intersectWithPlane( line, plane, _coeffs, angular_tolerance );
    }

    inline bool
    PlanePrimitive::intersectWithPlane( Eigen::Matrix<Scalar,-1,1>        & line
                                        , Eigen::Matrix<Scalar,4,1> const& plane_a
                                        , Eigen::Matrix<Scalar,4,1> const& plane_b
                                        , double                           angular_tolerance )
    {
        //planes shouldn't be parallel
        double test_cosine = plane_a.head<3>().dot(plane_b.head<3>());
        double upper_limit = 1 + angular_tolerance;
        double lower_limit = 1 - angular_tolerance;

        if ((test_cosine < upper_limit) && (test_cosine > lower_limit))
        {
            std::cerr<< "[" << __func__ << "] " << "Plane A and Plane B are Parallel" << std::endl;
            return (false);
        }

        if ((test_cosine > -upper_limit) && (test_cosine < -lower_limit))
        {
            std::cerr << "[" << __func__ << "] " << "Plane A and Plane B are Parallel" << std::endl;
            return (false);
        }

        Eigen::Matrix<Scalar,4,1> line_direction = plane_a.cross3(plane_b);
        line_direction.normalized();

        //construct system of equations using lagrange multipliers with one objective function and two constraints
        Eigen::Matrix<Scalar,-1,-1> langegrange_coefs(5,5);
        langegrange_coefs << 2,0,0,plane_a[0],plane_b[0],  0,2,0,plane_a[1],plane_b[1],  0,0,2, plane_a[2], plane_b[2], plane_a[0], plane_a[1] , plane_a[2], 0,0, plane_b[0], plane_b[1], plane_b[2], 0,0;

        Eigen::Matrix<Scalar,-1,1> b;
        b.resize(5);
        b << 0, 0, 0, -plane_a[3], -plane_b[3];

        //solve for the lagrange Multipliers
        Eigen::Matrix<Scalar,-1,1> x;
        x.resize(5);
        x = langegrange_coefs.colPivHouseholderQr().solve(b);

        line.resize(6);
        line.head<3>() = x.head<3>(); // the x[3] and x[4] are the values of the lagrange multipliers and are neglected
        line[3] = line_direction[0];
        line[4] = line_direction[1];
        line[5] = line_direction[2];

        return true;
    }
#endif
} //...ns GF2

#endif // __GF2_INC_PLANEPRIMITIVE_HPP__

#endif // __PLANEPRIMITIVE_H__

