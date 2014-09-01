#ifndef __GF2_LINEPRIMITIVE_H__
#define __GF2_LINEPRIMITIVE_H__

#include <fstream>

#include <Eigen/Dense>
#include "globfit2/primitives/primitive.h"

#ifdef GF2_USE_PCL
#   include "pcl/point_cloud.h"
#   include "pcl/sample_consensus/sac_model_line.h"
#   include "pcl/PointIndices.h"
#   include "pcl/visualization/pcl_visualizer.h"
#endif

namespace GF2
{
    //! \brief   Class to wrap a line with. Implements #pos() and #dir() functions, and a constructor taking a position and a direction.
    //!
    //!          Stores 3D position at the first three coeffs, and 3D direction at the second three.
    //! \warning Note, stores direction, NOT normal.
    class LinePrimitive : public ::GF2::Primitive<6>
    {
        public:
            // ____________________CONSTRUCT____________________
            //! \brief Inherited constructor from Primitive.
            using ::GF2::Primitive<6>::Primitive;

            //! \brief          Creates LinePrimitive from point on line and direction.
            //! \warning        NOT endpoints, use #fromEndPoints for that!
            //! \param[in] p0   Point on line.
            //! \param[in] dir  Line direction.
            LinePrimitive( Eigen::Matrix<Scalar,3,1> const& p0, Eigen::Matrix<Scalar,3,1> const& dir )
            {
                _coeffs.head   <3>( ) = p0;
                _coeffs.segment<3>(3) = dir.normalized();
            }

            static inline LinePrimitive
            fromEndPoints(  Eigen::Matrix<Scalar,3,1> p0, Eigen::Matrix<Scalar,3,1> p1 )
            {
                return LinePrimitive( p0, (p1-p0).normalized() );
            }

            // ____________________VIRTUALS____________________
            //! \brief  Compulsory virtual overload of position getter. The position of the line is the location stored at the first three coordinates of #_coeffs.
            //! \return The position of the line as a 3D Eigen::Vector.
            virtual Eigen::Matrix<Scalar,3,1> pos() const { return _coeffs.head   <3>( ); }
            //! \brief  Compulsory virtual overload of orientation getter. The orientation of the line is the direction stored at the second three coordinates of #_coeffs.
            //! \return The orientation of the line as a 3D Eigen::Vector.
            virtual Eigen::Matrix<Scalar,3,1> dir() const { return _coeffs.segment<3>(3); }

            //! \brief                  Snippet to get the normal of a line that lines on a plane specified by \p plane_normal.
            //! \tparam Scalar          Scalar type to use for calculations. Concept: float.
            //! \param[in] plane_normal Normal of plane, that needs to contain the normal of the line, that this function returns.
            //! \return                 Normal of the line lying in the plane specified by \p plane_normal.
            template <typename Scalar> inline Eigen::Matrix<Scalar,3,1>
            normal( Eigen::Matrix<Scalar,3,1> plane_normal = (Eigen::Matrix<Scalar,3,1>() << 0.,0.,1.).finished() ) const
            {
                Eigen::Matrix<Scalar,3,1> perp  = plane_normal * dir().dot( plane_normal );
                Eigen::Matrix<Scalar,3,1> par   = (dir() - perp).normalized();

                return par.cross( plane_normal ).normalized();
            } //...normal()

            //! \brief              Returns point to line distance.
            //! \param[in] point    Point to calculate distance from.
            //! \return             Distance from point to line.
            inline Scalar
            getDistance( Eigen::Matrix<Scalar,3,1> const& point ) const
            {
                return (this->pos() - point).cross( this->dir() ).norm();
            } //...getDistance()

            //! \brief                          Calculates the length of the line based on the points in \p cloud, masked by \p indices and the distance from point to line \p threshold.
            //!
            //!                                 The method calculates the inliers and selects the most far away point from #pos() in both directions.
            //! \tparam     _PointT             Type of endpoints returned in \p minMax. Concept: _PointContainerT::element_type::PointType == pcl::PointXYZRGB.
            //! \tparam     _PointContainerPtrT Type to store the points to select inliers from. Concept: pcl::PointCloud< _PointT >::Ptr.
            //! \param[out] minMax              Output container holding the two endpoints.
            //! \param[in]  cloud               Point container to look for inlier points in.
            //! \param[in]  threshold           A point in \p cloud is an inlier, if it is closer than this threshold. Usually the "scale".
            //! \param[in]  indices             Optional input to specify subset of points by indices.
            //! \return                         EXIT_SUCCESS
            //! \todo                           Detach from PCL.
            template <typename _PointT, typename _PointContainerPtrT>
            int
            getExtent( std::vector<_PointT>     & minMax
                     , _PointContainerPtrT        cloud
                     , const double               threshold = 0.01
                     , std::vector<int>    const* indices   = NULL ) const
            {
#           if GF2_USE_PCL
                pcl::PointIndices::Ptr                inliers( new pcl::PointIndices() );
                typename pcl::PointCloud<_PointT>::Ptr on_line_cloud( new pcl::PointCloud<_PointT>() );
                {
                    pcl::SampleConsensusModelLine<_PointT> sacline( cloud );
                    if ( indices )
                        sacline.setIndices( *indices );
                    sacline.selectWithinDistance(          coeffs(), threshold, inliers->indices      );
                    sacline.projectPoints       ( inliers->indices ,  coeffs(), *on_line_cloud, false );
                }

                if ( !on_line_cloud->size() )
                {
                    return EXIT_FAILURE;
                }

                float min_dist = 0.f, max_dist = 0.f;
                int   min_id   = 0  , max_id = 0;
                Eigen::Vector3f p0       = on_line_cloud->at(0).getVector3fMap();
                Eigen::Vector3f line_dir = this->dir(); line_dir.normalize();
                for ( size_t point_id = 1; point_id != on_line_cloud->size(); ++point_id )
                {
                    Eigen::Vector3f p1 = on_line_cloud->at(point_id).getVector3fMap();
                    Eigen::Vector3f p0p1 = p1-p0;
                    float dist = p0p1.dot( p0 + line_dir );
                    if ( dist < min_dist )
                    {
                        min_dist = dist;
                        min_id = point_id;
                    }
                    else if ( dist > max_dist )
                    {
                        max_dist = dist;
                        max_id = point_id;
                    }
                }

                // output
                minMax.resize( 2 );
                minMax[0] = on_line_cloud->at(min_id);
                minMax[1] = on_line_cloud->at(max_id);

                return EXIT_SUCCESS;
#           else
                std::cerr << "[" << __func__ << "]: " << "Needs PCL to work!!!!"
                return exit_FAILURE;
#           endif //...GF2_USE_PCL
            } //...draw()

            // ____________________DRAWING____________________
#if GF2_USE_PCL
            //! \brief              Draws a line using the two endpoints in \p ps on the visualizer \p vptr by the name \p name with colours specified by \p r, \p g, \p b on the viewport specified by \p viewport_id.
            //! \param[in] ps       Two endpoints.
            //! \param[in] vptr     Visualizer to draw on.
            static inline int
            draw( std::vector<pcl::PointXYZ>              const& ps
                , pcl::visualization::PCLVisualizer::Ptr         vptr
                , std::string                             const  name
                , double r, double g, double b
                , int                                     const  viewport_id = 0 )
            {
                if ( ps.size() != 2 )
                {
                    std::cerr << "[" << __func__ << "]: " << "ps.size() " << ps.size() << " != 2 ...exiting\n"; fflush(stderr);
                    return EXIT_FAILURE;
                }
                return vptr->addLine( ps[0]
                                    , ps[1]
                                    , r, g, b
                                    , name
                                    , viewport_id );
            } //...draw()

            //! \brief              Draws a line as long as it finds 2 points inside its radius. The radius is doubled 10 times until it gives up looking for more points.
            //! \tparam PointsPtrT  Holds the points. Is usually a pcl::PointCloud::Ptr of some sort. Assumed to be a pointer.
            //! \tparam PointT      Type of point that is stored in the PointCloud.
            //! \param  stretch     Multiplies the final length with this number, so that the line is a bit longer than the distance between its extrema.
            template <class PointsPtrT, class PointT = typename PointsPtrT::element_type::PointType, class Scalar = float>
            static inline int
            draw( LinePrimitive                    const& line
                , PointsPtrT                              cloud
                , Scalar                                  radius
                , std::vector<int>                 const* indices
                , pcl::visualization::PCLVisualizer::Ptr  v
                , std::string                             plane_name
                , double                                  r, double g, double b
                , int                                     viewport_id = 0
                , Scalar                                  stretch     = 1.f
                )
            {
                int err     = EXIT_SUCCESS;

                std::vector<PointT> minMax;
                int     it          = 0;
                int     max_it      = 10;
                Scalar  tmp_radius  = radius;
                do
                {
                    err = line.getExtent( minMax
                                        , cloud
                                        , tmp_radius
                                        , indices   );
                    tmp_radius *= 2.f;
                } while ( (minMax.size() < 2) && (++it < max_it) );

                if ( it >= max_it )
                {
                    std::cerr << "[" << __func__ << "]: " << "line.getExtent exceeded max radius increase iteration count...not drawing" << line.toString() << std::endl;
                    return err;
                }

                std::vector<pcl::PointXYZ> ps;

                Eigen::Map<Eigen::Vector3f> p0 = minMax[0].getVector3fMap();
                Eigen::Map<Eigen::Vector3f> p1 = minMax[1].getVector3fMap();
                Eigen::Matrix<Scalar,3,1> diff = p1 - p0;
                Scalar half_stretch = Scalar(1) + (stretch-Scalar(1)) / Scalar(2.);
                Eigen::Matrix<Scalar,3,1> p1_final  = p0 + diff * half_stretch;
                Eigen::Matrix<Scalar,3,1> p0_final  = p1 - diff * half_stretch;

                pcl::PointXYZ pnt;
                pnt.x = p0_final(0); pnt.y = p0_final(1); pnt.z = p0_final(2);
                ps.push_back( pnt );
                pnt.x = p1_final(0); pnt.y = p1_final(1); pnt.z = p1_final(2);
                ps.push_back( pnt );

                err += draw( ps, v, plane_name, r, g, b, viewport_id );

                return err;
            }
#endif // GF2_USE_PCL

            // ____________________DEPRECATED____________________
            //! \deprecated Decides, if two lines are different up to a position and and an angle threshold.
            static bool
            different( LinePrimitive const& me, LinePrimitive const& other, Scalar pos_diff, Scalar ang_diff )
            {
                // pos
                if ( (me.pos() - other.pos()).norm() > pos_diff ) return true;

                // angle
                Scalar angle = acos( me.dir().dot(other.dir()) ); if ( angle != angle ) angle = static_cast<Scalar>(0);
                return ( angle > ang_diff );
            } //...different()

#if 0
            static Scalar
            distanceToPoint( Eigen::Matrix<Scalar,Dim,1> line, Eigen::Matrix<Scalar,4,1> point, bool squared = false )
            {
                  // Obtain the line point and direction
                  Eigen::Matrix<Scalar,4,1> line_pt ( line[0], line[1], line[2], 0 );
                  Eigen::Matrix<Scalar,4,1> line_dir( line[3], line[4], line[5], 0 );
                  line_dir.normalize ();

                  // Calculate the distance from the point to the line
                  // D = ||(P2-P1) x (P1-P0)|| / ||P2-P1|| = norm (cross (p2-p1, p2-p0)) / norm(p2-p1)
                  // Need to estimate sqrt here to keep MSAC and friends general
                  if ( squared )    return (line_pt - point).cross3(line_dir).squaredNorm();
                  else              return (line_pt - point).cross3(line_dir).norm();
            } //...distanceToPoint()

            inline Scalar distanceToPoint( Eigen::Matrix<Scalar,4,1> const& point, bool squared = false ) const { return distanceToPoint( _coeffs,   point, squared ); }

            static inline Scalar
            point3Distance( Eigen::Matrix<Scalar,3,1> const& point, Eigen::Matrix<Scalar,-1,1> const& line )
            {
                return (line.template head<3>() - point).cross( line.template segment<3>(3) ).norm();
            }

            static inline Scalar
            point4Distance( Eigen::Matrix<Scalar,4,1> const& point, Eigen::Matrix<Scalar,-1,1> const& line )
            {
                  // Calculate the distance from the point to the line
                  // D = ||(P2-P1) x (P1-P0)|| / ||P2-P1|| = norm (cross (p2-p1, p2-p0)) / norm(p2-p1)
                  return (line.template head<3>() - point.head<3>()).cross(line.template segment<3>(3)).norm();
            }

            // non-static relays
            inline Scalar point3Distance ( Eigen::Matrix<Scalar,3,1> const& point                       ) const { return point3Distance (   point, _coeffs          ); }
            inline Scalar point4Distance ( Eigen::Matrix<Scalar,4,1> const& point                       ) const { return point4Distance (   point, _coeffs          ); }
#endif


    }; // LinePrimitive

} // ns GF2

#endif // __GF2_LINEPRIMITIVE_H__
