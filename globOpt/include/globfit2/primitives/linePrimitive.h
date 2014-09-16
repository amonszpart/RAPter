#ifndef __GF2_LINEPRIMITIVE_H__
#define __GF2_LINEPRIMITIVE_H__

#include <fstream>
#include <memory> // shared_ptr

#include <Eigen/Dense>
#include "globfit2/primitives/primitive.h" // inherit
#include "globfit2/primitives/taggable.h"  // inherit

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
    class LinePrimitive : public ::GF2::Primitive<2,6>, public ::GF2::Taggable
    {
            typedef ::GF2::Primitive<2,6> ParentT;
        public:
            //! \brief Defines the tags (ids) that this primitive can manage using setTag and getTag functions.
            enum TAGS {
                GID        = 0  //!< group id             - which group this primitive is supposed to explain
                , DIR_GID  = 1  //!< direction group id   - which group this primitive got it's direction from
                , CHOSEN   = 2  //!< an additional flag to store, if this is part of a solution.
//                , USER_ID1 = 10 //!< additional flag to store processing attributes (values only in the generation scope)
//                , USER_ID2 = 11 //!< additional flag to store processing attributes (values only in the generation scope)
//                , USER_ID3 = 12 //!< additional flag to store processing attributes (values only in the generation scope)
//                , USER_ID4 = 13 //!< additional flag to store processing attributes (values only in the generation scope)
//                , USER_ID5 = 14 //!< additional flag to store processing attributes (values only in the generation scope)
            }; //...TAGS

            typedef ParentT::Scalar Scalar;

            // ____________________CONSTRUCT____________________
            LinePrimitive() : ParentT() {}

            //! \brief Constructor that takes raw data in Eigen format as input.
            LinePrimitive( Eigen::Matrix<Scalar,Dim,1> coeffs ) : ParentT( coeffs ) {}

            //! \brief Constructor that takes raw data in std::vector format as input.
            LinePrimitive( std::vector<Scalar> const& coeffs ) : ParentT( coeffs ) {}

            //! \brief          Creates LinePrimitive from point on line and direction.
            //! \warning        NOT endpoints, use #fromEndPoints for that!
            //! \param[in] p0   Point on line.
            //! \param[in] dir  Line direction.
            LinePrimitive( Eigen::Matrix<Scalar,3,1> const& p0, Eigen::Matrix<Scalar,3,1> const& dir )
            {
                _coeffs.head   <3>( ) = p0;
                _coeffs.segment<3>(3) = dir.normalized();
            }

            LinePrimitive( Eigen::Matrix<Scalar,3,1> const& centroid, Eigen::Matrix<Scalar,3,1> const& eigen_values, Eigen::Matrix<Scalar, 3, 3> const& eigen_vectors )
            {
                // get eigen vector for biggest eigen value
                const int max_eig_val_id = std::distance( eigen_values.data(), std::max_element( eigen_values.data(), eigen_values.data()+3 ) );
                // set position
                _coeffs.template head<3>() = centroid;
                // set direction
                _coeffs.template segment<3>(3) = eigen_vectors.col(max_eig_val_id).normalized();
            }

            static inline LinePrimitive
            fromEndPoints(  Eigen::Matrix<Scalar,3,1> p0, Eigen::Matrix<Scalar,3,1> p1 )
            {
                return LinePrimitive( p0, (p1-p0).normalized() );
            }

            /*! \brief Called from \ref GF2::CandidateGenerator::generate to create a new candidate from this and \p other.
             *         Create back-rotated version of the other at this position.
             *  \param[out] out              Created output primitive.
             *  \param[in]  other            Primitive, who's direction we want to use to generate something at our location.
             *  \param[in]  closest_angle_id The id of closest perfect angle between us and \p other. This we want to rotate back by.
             *  \param[in]  angles           List of desired angles, one of its entries is referenced by closest_angle_id.
             *  \param[in]  angle_multiplier A coefficient to multiply the angle by. Set to -1, if we want to rotate back.
             *  \return \p out is only valid, if return true. We might decide, that there's nothing to create, in which case we return false.
             */
            template <class _AngleContainerT>
            inline bool generateFrom( LinePrimitive         & out
                                    , LinePrimitive    const& other
                                    , int              const  closest_angle_id
                                    , _AngleContainerT const& angles
                                    , Scalar           const  angle_multiplier = Scalar(1.)
                                    ) const
            {
                Scalar const angle = angles[ closest_angle_id ] * angle_multiplier;

                out = LinePrimitive( /*  position: */ this->pos()
                                   , /* direction: */ Eigen::AngleAxisf(angle, Eigen::Matrix<Scalar,3,1>::UnitZ())
                                                      * other.dir()
                                   );
                // copy position id from self
                out.setTag( GID, this->getTag(GID) );
                // copy direction id from the other
                out.setTag( DIR_GID, other.getTag(DIR_GID) );
                // erase chosen tag - this is a new candidate
                out.setTag( CHOSEN , -1 );

                return true;
            } //...generateFrom

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
            inline Eigen::Matrix<Scalar,3,1>
            normal( Eigen::Matrix<Scalar,3,1> plane_normal = (Eigen::Matrix<Scalar,3,1>() << 0.,0.,1.).finished() ) const
            {
                Eigen::Matrix<Scalar,3,1> perp  = plane_normal * dir().dot( plane_normal );
                Eigen::Matrix<Scalar,3,1> par   = (dir() - perp).normalized();

                return par.cross( plane_normal ).normalized();
            } //...normal()

            // _______________________IO_______________________
            /*! \brief Used in \ref io::readPrimitives to determine how many floats to parse from one entry.
             */
            static inline int getFileEntryLength() { return 6; }

            /*! \brief Output <x0,n>, location and normal for the plane to a string that does *not* have an endline at the end.
             */
            inline std::string toFileEntry() const
            {
                Eigen::Matrix<Scalar,3,1> pos = this->pos();
                Eigen::Matrix<Scalar,3,1> nrm = this->normal();
                char line[1024];
                sprintf( line, "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,"
                         , pos(0), pos(1), pos(2)
                         , nrm(0), nrm(1), nrm(2) );
                return std::string( line );
            }

            /*! \brief              Used in \ref io::readPrimitives to create the primitive from the read floats.
             *  \param[in] entries  Contains <x0,n> to create the plane from.
             */
            static inline LinePrimitive fromFileEntry( std::vector<Scalar> const& entries )
            {
                return LinePrimitive( Eigen::Map<const Eigen::Matrix<Scalar,3,1> >( entries.data()  , 3 ),
                                       Eigen::Map<const Eigen::Matrix<Scalar,3,1> >( entries.data()+3, 3 ).cross(Eigen::Matrix<Scalar,3,1>::UnitZ()) );
            } //...fromFileEntry()

            // ____________________GEOMETRY____________________

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
            //! \tparam     _PointT             Type of point stored in \p cloud. Concept: \ref GF2::PointPrimitve.
            //! \tparam     _Scalar             Precision of data in \p minMax.
            //! \tparam     _PointContainerPtrT Type to store the points to select inliers from. Concept: "std::vector<PointPrimitive>*" .
            //! \param[out] minMax              Output container holding the two endpoints.
            //! \param[in]  cloud               Point container to look for inlier points in.
            //! \param[in]  threshold           A point in \p cloud is an inlier, if it is closer than this threshold. Usually the "scale".
            //! \param[in]  indices             Optional input to specify subset of points by indices.
            //! \return                         EXIT_SUCCESS
            template <typename _PointT, typename _PointContainerT>
            int
            getExtent( std::vector<Eigen::Matrix<Scalar,3,1> >      & minMax
                     , _PointContainerT                         const& cloud
                     , double                                   const  threshold   = 0.01
                     , std::vector<int>                         const* indices_arg = NULL ) const
            {
                typedef Eigen::Matrix<Scalar,3,1> Position;

                // select inliers
                std::vector<int> inliers;
                inliers.reserve( cloud.size() );
                const int stop_at = indices_arg ? indices_arg->size() : cloud.size();
                for ( int i = 0; i != stop_at; ++i )
                {
                    const int pid = indices_arg ? (*indices_arg)[i] : i;
                    if ( this->getDistance( cloud[pid].template pos() ) < threshold )
                        inliers.push_back( pid );
                }

                // check size
                if ( !inliers.size() ) return EXIT_FAILURE;

                // project cloud
                std::vector<Position> on_line_cloud;
                for ( int pid_id = 0; pid_id != inliers.size(); ++pid_id )
                {
                    const int pid = inliers[ pid_id ];
                    on_line_cloud.push_back( this->projectPoint(cloud[pid].template pos()) );
                }


                // select min-max inliier from projected cloud
                Scalar min_dist = 0.f, max_dist = 0.f;
                int     min_id   = 0  , max_id   = 0;
                Position    p0       = on_line_cloud[0];
                Position    line_dir = this->dir();

                for ( size_t point_id = 1; point_id != on_line_cloud.size(); ++point_id )
                {
                    Position const& p1   = on_line_cloud[ point_id ];
                    Position        p0p1 = p1-p0;
                    float           dist = p0p1.dot( p0 + line_dir );
                    if ( dist < min_dist )
                    {
                        min_dist = dist;
                        min_id   = point_id;
                    }
                    else if ( dist > max_dist )
                    {
                        max_dist = dist;
                        max_id   = point_id;
                    }
                }

                // output
                minMax.resize( 2 );
                minMax[0] = on_line_cloud[ min_id ];
                minMax[1] = on_line_cloud[ max_id ];

                return EXIT_SUCCESS;
            } //...draw()

            inline Eigen::Matrix<Scalar,3,1>
            projectPoint( Eigen::Matrix<Scalar,3,1> const& point ) const
            {
                Eigen::Matrix<Scalar,4,1> line_pt ; line_pt  << this->pos(), 0;
                Eigen::Matrix<Scalar,4,1> line_dir; line_dir << this->dir(), 0;
                Eigen::Matrix<Scalar,4,1> pt      ; pt       << point, 0;

                Scalar k = (pt.dot (line_dir) - line_pt.dot (line_dir)) / line_dir.dot (line_dir);

                Eigen::Matrix<Scalar,4,1> pp = line_pt + k * line_dir;
                return pp.template head<3>();
            } // projectPoint

            // ____________________DRAWING____________________
#if GF2_USE_PCL
            /*! \brief              Draws a line using the two endpoints in \p ps on the visualizer \p vptr by the name \p name with colours specified by \p r, \p g, \p b on the viewport specified by \p viewport_id.
             *  \param[in] ps       Two endpoints.
             *  \param[in] vptr     Visualizer to draw on.
             */
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

            /*! \brief              Draws a line as long as it finds 2 points inside its radius. The radius is doubled 10 times until it gives up looking for more points.
             *  \tparam _PointPrimitiveT   Concept: \ref GF2::PointPrimitive.
             *  \tparam _PointContainerT   Concept: std::vector< _PointPrimitiveT >.
             *  \tparam _IndicesContainerT Concept: std::vector<int>.
             *  \param  stretch     Multiplies the final length with this number, so that the line is a bit longer than the distance between its extrema.
             */
            template <class _PointPrimitiveT, class _PointContainerT>
            static inline int
            draw( LinePrimitive                    const& line
                , _PointContainerT                 const& cloud
                , Scalar                           const  radius
                , std::vector<int>                 const* indices
                , pcl::visualization::PCLVisualizer::Ptr  v
                , std::string                      const  plane_name
                , double                           const  r
                , double                           const  g
                , double                           const  b
                , int                              const  viewport_id = 0
                , Scalar                          const  stretch     = Scalar( 1. )
                )
            {
                typedef Eigen::Matrix<Scalar,3,1> Position;

                int err     = EXIT_SUCCESS;

                std::vector< Position > minMax;
                int      it          = 0;
                int      max_it      = 10;
                Scalar  tmp_radius  = radius;
                do
                {
                    err = line.getExtent<_PointPrimitiveT>( minMax
                                                          , cloud
                                                          , tmp_radius
                                                          , indices   );
                    tmp_radius *= 2.f;
                } while ( (minMax.size() < 2) && (++it < max_it) );

                if ( it >= max_it )
                {
                    //std::cerr << "[" << __func__ << "]: " << "line.getExtent exceeded max radius increase iteration count...not drawing " << line.toString() << std::endl;
                    std::cerr << "[" << __func__ << "]: " << "line.getExtent exceeded max radius increase iteration count...drawing unit " << line.toString() << std::endl;
                    minMax.resize(2);
                    minMax[0] = line.pos();
                    minMax[1] = line.pos() + line.dir() / 10.;
                    //return err;
                }

                std::vector<pcl::PointXYZ> ps;

                Position const& p0 = minMax[0];
                Position const& p1 = minMax[1];
                Position diff = p1 - p0;
                Scalar half_stretch = Scalar(1) + (stretch-Scalar(1)) / Scalar(2.);
                Position p1_final  = p0 + diff * half_stretch;
                Position p0_final  = p1 - diff * half_stretch;

                pcl::PointXYZ pnt;
                pnt.x = p0_final(0); pnt.y = p0_final(1); pnt.z = p0_final(2);
                ps.push_back( pnt );
                pnt.x = p1_final(0); pnt.y = p1_final(1); pnt.z = p1_final(2);
                ps.push_back( pnt );

                err += draw( ps, v, plane_name, r, g, b, viewport_id );

                return err;
            }
#endif // GF2_USE_PCL

    }; // LinePrimitive

} // ns GF2

#endif // __GF2_LINEPRIMITIVE_H__
