#ifndef __GF2_LINEPRIMITIVE_H__
#define __GF2_LINEPRIMITIVE_H__

#include <fstream>
#include <memory> // shared_ptr

#include <Eigen/Dense>
#include "globfit2/optimization/energyFunctors.h"   // pointToFiniteLineDistanceFunctor
#include "globfit2/primitives/primitive.h"          // inherit
#include "globfit2/processing/util.hpp"             // getPopulationOf()

#ifdef GF2_USE_PCL
#   include "pcl/point_cloud.h"
#   include "pcl/sample_consensus/sac_model_line.h"
#   include "pcl/PointIndices.h"
#   include "pcl/visualization/pcl_visualizer.h"
#   include "pcl/surface/concave_hull.h"
#endif

namespace GF2
{
    //! \brief   Class to wrap a line with. Implements #pos() and #dir() functions, and a constructor taking a position and a direction.
    //!
    //!          Stores 3D position at the first three coeffs, and 3D direction at the second three.
    //! \warning Note, stores direction, NOT normal.
    class LinePrimitive : public ::GF2::Primitive<2,6>
    {
            typedef ::GF2::Primitive<2,6>   ParentT;
        public:
            typedef ParentT::Scalar         Scalar;
            typedef ParentT::Position       Position;
            typedef ParentT::ExtremaT       ExtremaT;

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
            LinePrimitive( Position const& p0, Direction const& dir )
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
                                    , bool             const  verbose = false
                                    ) const
            {
                typedef Eigen::Matrix<Scalar,3,1> Direction;

                Scalar const angle = angles[ closest_angle_id ] * angle_multiplier;
                Direction d0 = Eigen::AngleAxisf( angle, Direction::UnitZ()) * other.dir(),
                          d1 = Eigen::AngleAxisf(-angle, Direction::UnitZ()) * other.dir();

                bool revert = false;
                if ( angleInRad(this->dir(),d0) > angleInRad(this->dir(),d1) )
                     revert = true;

                out = LinePrimitive( /*  position: */ this->pos()
                                   , /* direction: */ revert ? d1 : d0
                                   );

                out.copyTagsFrom( *this ); // added by Aron on 3/1/2015
                // copy position id from self
                //out.setTag( GID, this->getTag(GID) );
                // copy direction id from the other
                out.setTag( TAGS::DIR_GID, other.getTag(TAGS::DIR_GID) );
                // erase chosen tag - this is a new candidate
                out.setTag( TAGS::STATUS, STATUS_VALUES::UNSET );

                if ( verbose ) std::cout << "[" << __func__ << "]: " << "angle= " << angle << ", closest_angle_id: " << closest_angle_id << ", mult:  " << angle_multiplier << std::endl;

                return true;
            } //...generateFrom

            template <typename DerivedT> static
            int generateFrom( LinePrimitive      & out
                            , DerivedT      const& normal
                            , Scalar        const  distanceFromOrigin );

            // ____________________VIRTUALS____________________
            //! \brief  Compulsory virtual overload of position getter. The position of the line is the location stored at the first three coordinates of #_coeffs.
            //! \return The position of the line as a 3D Eigen::Vector.
            //virtual Eigen::Matrix<Scalar,3,1> pos() const { return _coeffs.head   <3>( ); }
            inline typename Eigen::Matrix<Scalar,Dim,1>::ConstFixedSegmentReturnType<3>::Type pos() const { return _coeffs.template head<3>(); }
            //! \brief  Compulsory virtual overload of orientation getter. The orientation of the line is the direction stored at the second three coordinates of #_coeffs.
            //! \return The orientation of the line as a 3D Eigen::Vector.
            //virtual Eigen::Matrix<Scalar,3,1> dir() const { return _coeffs.segment<3>(3); }
            inline typename Eigen::Matrix<Scalar,Dim,1>::ConstFixedSegmentReturnType<3>::Type dir() const { return _coeffs.template segment<3>(3); }

            //! \brief                  Snippet to get the normal of a line that lines on a plane specified by \p plane_normal.
            //! \tparam Scalar          Scalar type to use for calculations. Concept: float.
            //! \param[in] plane_normal Normal of plane, that needs to contain the normal of the line, that this function returns.
            //! \return                 Normal of the line lying in the plane specified by \p plane_normal.
            inline Direction
            normal( Direction plane_normal = (Direction() << 0.,0.,1.).finished() ) const
            {
                Direction perp  = plane_normal * dir().dot( plane_normal );
                Direction par   = (dir() - perp).normalized();

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

            /*! \brief              Helps calculating line to line distance.
             *  \param[in] extrema  Extrema of this primitive.
             *  \param[in] pnt      One of the extrema of the other primitive.
             *  \return             Distance from point to plane.
             */
            inline Scalar
            getFiniteDistance( ExtremaT const& extrema, Position const& pnt ) const
            {
                return MyPointFiniteLineDistanceFunctor::eval( extrema, *this, pnt );
            }

            /*! \brief                          Calculates the length of the line based on the points in \p cloud, masked by \p indices and the distance from point to line \p threshold.
             *
             *                                  The method calculates the inliers and selects the most far away point from #pos() in both directions.
             *  \tparam     _PointT             Type of point stored in \p cloud. Concept: \ref GF2::PointPrimitve.
             *  \tparam     _Scalar             Precision of data in \p minMax.
             *  \tparam     _PointContainerPtrT Type to store the points to select inliers from. Concept: "std::vector<PointPrimitive>*" .
             *  \param[out] minMax              Output container holding the two endpoints.
             *  \param[in]  cloud               Point container to look for inlier points in.
             *  \param[in]  threshold           A point in \p cloud is an inlier, if it is closer than this threshold. Usually the "scale".
             *  \param[in]  indices             Optional input to specify subset of points by indices.
             *  \return                         EXIT_SUCCESS
             */
            template <typename _PointT, class _IndicesContainerT, typename _PointContainerT>
            int
            getExtent( ExtremaT                                      & minMax
                     , _PointContainerT                         const& cloud
                     , double                                   const  threshold            = 0.01
                     , _IndicesContainerT                       const* indices_arg          = NULL
                     , bool                                     const  force_axis_aligned   = false ) const;
#if 0
            {
                typedef Eigen::Matrix<Scalar,3,1> Position;

                if ( this->_extents.isUpdated() )
                {
                    minMax = _extents.get();
                    return 0;
                }

                // select inliers
                std::vector<PidT> inliers;
                inliers.reserve( cloud.size() );
                const PidT stop_at = indices_arg ? indices_arg->size() : cloud.size();
                for ( PidT i = 0; i != stop_at; ++i )
                {
                    const PidT pid = indices_arg ? (*indices_arg)[i] : i;
                    if ( this->getDistance( cloud[pid].template pos() ) < threshold )
                        inliers.push_back( pid );
                }

                // check size
                if ( !inliers.size() )
                {
                    std::cerr << "no inliers for primitive gid"
                              << this->getTag( TAGS::GID ) << ", did "
                              << this->getTag( TAGS::DIR_GID ) << std::endl;
                              return EXIT_FAILURE;
                }

                // project cloud
                std::vector<Position> on_line_cloud;
                for ( PidT pid_id = 0; pid_id != inliers.size(); ++pid_id )
                {
                    const PidT pid = inliers[ pid_id ];
                    on_line_cloud.push_back( this->projectPoint(cloud[pid].template pos()) );
                }


                // select min-max inliier from projected cloud
                Scalar min_dist = 0.f, max_dist = 0.f;
                PidT     min_id   = 0  , max_id   = 0;
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

                this->_extents.update( minMax );

                return EXIT_SUCCESS;
            } //...getExtent()
#endif
    public:

            /*! \brief Calculates size, a bit smarter, than taking the area of #getExtent().
             *  \tparam MatrixDerived   Concept: Eigen::Matrix<_Scalar,-1,1>.
             */
            template <class _IndicesContainerT, typename MatrixDerived, typename _Scalar, class _PointContainerT > inline
            MatrixDerived& getSpatialSignificance( MatrixDerived& in, _PointContainerT const& points, _Scalar const /*scale*/
                                                 , _IndicesContainerT *indices = NULL, bool return_squared = false ) const
            {
                _IndicesContainerT tmp_population,
                                   *pop           = &tmp_population;
                if ( !indices )
                    processing::getPopulationOf( tmp_population, this->getTag(TAGS::GID), points );
                else
                    pop = indices;

                if ( !(pop->size()) )
                {
                    std::cerr << "[" << __func__ << "]: " << "_____________NO points in primitive!!!!_____________" << std::endl;
                    in.setConstant( _Scalar(-1.) );
                }
                in.setConstant( _Scalar(0.) );

#if 1 // biggest eigen value
                Eigen::Matrix<_Scalar,3,1> eigen_values;
                Eigen::Matrix<_Scalar,3,3> eigen_vectors;
                processing::eigenDecomposition( eigen_values, eigen_vectors, points, pop );
                if ( return_squared )
                    in(0) = eigen_values( 0 );
                else
                    in(0) = std::sqrt( eigen_values(0) );
#else // variance
                Eigen::Matrix<_Scalar,3,1> centroid = processing::getCentroid<_Scalar>( points, &population );
                for ( size_t pid_id = 0; pid_id != population.size(); ++pid_id )
                {
                    const int pid = population[pid_id];
                    std::cout << "[" << __func__ << "]: " << "\tadding " \
                    "points[" << pid <<"].pos() (" << points[pid].template pos().transpose()
                              << " - " << centroid.transpose() << ").squaredNorm(): "
                              << (points[pid].template pos() - centroid).squaredNorm() << "\n";
                    in(0) += (points[pid].template pos() - centroid).squaredNorm();
                    std::cout << "\tin0 is now " << in(0) << std::endl;
                }
                in(0) /= _Scalar( population.size() );
#endif

                return in;
            } //...getSpatialSignificance()

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
            template <class _PointPrimitiveT, class _IndicesContainerT, class _PointContainerT>
            static int
            draw( LinePrimitive                    const& line
                , _PointContainerT                 const& cloud
                , Scalar                           const  radius
                , _IndicesContainerT               const* indices
                , pcl::visualization::PCLVisualizer::Ptr  v
                , std::string                      const  plane_name
                , double                           const  r
                , double                           const  g
                , double                           const  b
                , int                              const  viewport_id = 0
                , Scalar                           const  stretch     = Scalar( 1. )
                , int                              const  draw_mode   = 0               // this is needed in plane
                , float                            const  hull_alpha  = 0. // this is needed in plane
                );

            /*! \brief Extract convec hull and display
             *  \tparam PointsT Concept: std::vector<pcl::PointXYZ>.
             */
            template < class PclCloudT, class _PointContainerT, class _IndicesContainerT> static inline int
            getHull( PclCloudT                & plane_polygon_cloud // pcl::PointCloud<pcl::PointXYZ>::Ptr plane_polygon_cloud_ptr( new pcl::PointCloud<pcl::PointXYZ> );
                   , LinePrimitive       const& plane
                   , _PointContainerT    const& points
                   , _IndicesContainerT  const* indices
                   , float               const  alpha = 2.f
                   , pcl::PolygonMesh         * out_mesh = NULL
                   )
            {
                typedef typename PclCloudT::PointType PclPointT;

                pcl::ConcaveHull<PclPointT>              concave_hull;                                   // object
                //typename pcl::PointCloud<PclPointT>::Ptr cloud_hull( new pcl::PointCloud<PclPointT>() );
                typename pcl::PointCloud<PclPointT> cloud_hull;
                typename pcl::PointCloud<PclPointT> cloud_projected;
                std::vector<pcl::Vertices>                   polygons;                               // output list indexing the points from cloud_hull, in 2D this is size 1
                //pcl::PointCloud<PclPointT>::Ptr plane_polygon_cloud_ptr( new pcl::PointCloud<PclPointT> );

                cloud_projected.resize(indices->size());

                // get assigned points, project them to the plane and store as PCL cloud
                PidT i = 0;
                for ( typename _IndicesContainerT::const_iterator it = indices->begin(); it != indices->end(); ++i, ++it )
                {
                    cloud_projected.at(i).getVector3fMap() = plane.projectPoint(points[*it].pos()).template cast<float>();
                }

                concave_hull.setAlpha( alpha );
                concave_hull.setInputCloud( cloud_projected.makeShared() );
                concave_hull.reconstruct( cloud_hull, polygons );
                PidT max_size = 0, max_id = 0;
                for ( i = 0; i != polygons.size(); ++i )
                {
                    if ( polygons[i].vertices.size() >max_size)
                    {
                        max_size = polygons[i].vertices.size();
                        max_id = i;
                    }
                }

                plane_polygon_cloud.resize( polygons[max_id].vertices.size() );
                i = 0;
                for ( std::vector<uint32_t>::const_iterator it  = polygons[max_id].vertices.begin();
                                                            it != polygons[max_id].vertices.end(); ++i, ++it )
                {
                    plane_polygon_cloud.at( i ) = cloud_hull.at(*it);
                }

                if ( out_mesh )
                {
                    // Perform reconstruction
                    out_mesh->polygons = polygons;

                    // Convert the PointCloud into a PCLPointCloud2
                    pcl::toPCLPointCloud2 (cloud_hull, out_mesh->cloud);
                }


                return plane_polygon_cloud.size();
            }
#endif // GF2_USE_PCL

            inline bool gidUnset() const { return this->getTag(LinePrimitive::TAGS::GID) == LONG_VALUES::UNSET; }

    }; // LinePrimitive

} // ns GF2

#endif // __GF2_LINEPRIMITIVE_H__
