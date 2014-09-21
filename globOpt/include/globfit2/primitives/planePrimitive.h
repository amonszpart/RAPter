#ifndef __GF2_PLANEPRIMITIVE_H__
#define __GF2_PLANEPRIMITIVE_H__

#include <Eigen/Dense>
#include "globfit2/primitives/primitive.h"
//#include "pcltools/util.hpp"
#include "globfit2/processing/util.hpp" // pca

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
    class PlanePrimitive : public ::GF2::Primitive<3,6>
    {
            typedef ::GF2::Primitive<3,6> ParentT;
        public:

            typedef ParentT::Scalar Scalar;

            //! \brief Inherited constructor from Primitive.
            // using ::GF2::Primitive<4>::Primitive;

            // ____________________CONSTRUCT____________________

            PlanePrimitive() : ParentT() {}

            //! \brief Constructor that takes raw data in Eigen format as input.
            PlanePrimitive( Eigen::Matrix<Scalar,Dim,1> coeffs ) : ParentT( coeffs ) {}

            //! \brief Constructor that takes raw data in std::vector format as input.
            PlanePrimitive( std::vector<Scalar> const& coeffs ) : ParentT( coeffs ) {}

            //! \brief           Creates PlanePrimitive from point on plane and direction.
            //! \param[in] p0    Point on plane.
            //! \param[in] dir   Plane normal.
            PlanePrimitive( Eigen::Matrix<Scalar,3,1> pnt, Eigen::Matrix<Scalar,3,1> normal );

            /*! \brief Creates plane primitive from an ordered local frame
             */
            PlanePrimitive( Eigen::Matrix<Scalar,3,1> const& centroid, Eigen::Matrix<Scalar,3,1> const& eigen_values, Eigen::Matrix<Scalar, 3, 3> const& eigen_vectors )
            {
                // get eigen vector for biggest eigen value
                const int min_eig_val_id = std::distance( eigen_values.data(), std::min_element( eigen_values.data(), eigen_values.data()+3 ) );
                // set position
                _coeffs.template head<3>() = centroid;
                // set direction
                _coeffs.template segment<3>(3) = eigen_vectors.col(min_eig_val_id).normalized();
            }

            /*! \brief Called from \ref GF2::CandidateGenerator::generate to create a new candidate from this and \p other.
             *         Create back-rotated version of the other at this position.
             *  \param[out] out              Created output primitive.
             *  \param[in]  other            Primitive, who's direction we want to use to generate something at our location.
             *  \param[in]  closest_angle_id The id of closest perfect angle between us and \p other. This we want to rotate back by.
             *  \param[in]  angles           List of desired angles, one of its entries is referenced by closest_angle_id.
             *  \return \p out is only valid, if return true. We might decide, that there's nothing to create, in which case we return false.
             */
            template <class _AngleContainerT>
            inline bool generateFrom( PlanePrimitive         & out
                                    , PlanePrimitive    const& other
                                    , int               const  closest_angle_id
                                    , _AngleContainerT  const& angles
                                    , Scalar            const  angle_multiplier = Scalar(1.)
                                    ) const
            {
                // if not 0 or M_PI, meaning not parallel
                if ( (closest_angle_id != 0) && (closest_angle_id != angles.size()-1) )
                {
                    //std::cerr << "[" << __func__ << "]: " << "rotating plane by angle " << angles[closest_angle_id] << std::endl;
                    //return false;
                    //Scalar const angle = angles[ closest_angle_id ];
                    Scalar angle = angles[ closest_angle_id ];
                    out = PlanePrimitive( /*  position: */ this->pos()
                                        , /* direction: */ Eigen::AngleAxisf( angle, other.dir().cross(dir()) ) * other.dir()
                                        );
                }
                else
                {
                    //Scalar const angle = angles[ closest_angle_id ];
                    out = PlanePrimitive( /*  position: */ this->pos()
                                        , /* direction: */ other.dir()
                                        );
                }
                // copy position id from self
                out.setTag( GID, this->getTag(GID) );
                // copy direction id from the other
                out.setTag( DIR_GID, other.getTag(DIR_GID) );
                // erase chosen tag - this is a new candidate
                out.setTag( STATUS, TAG_UNSET );

                return true;
            } //...generateFrom

            // ____________________VIRTUALS____________________
            //! \brief  Compulsory virtual overload of position getter. The position of the plane is calculated on the fly from the formula N . x0 + d = 0.
            //! \return The position of the plane as a 3D Eigen::Vector.
            //virtual Eigen::Matrix<Scalar,3,1> pos() const { return _coeffs.head<3>() * (-1.f*_coeffs(3)); }
            virtual Eigen::Matrix<Scalar,3,1> pos() const { return _coeffs.head<3>(); }
            //! \brief  Compulsory virtual overload of orientation getter. The orientation of the plane is the normal stored at the first three coordinates of #_coeffs.
            //! \return The normal of the plane as a 3D Eigen::Vector.
            //virtual Eigen::Matrix<Scalar,3,1> dir() const { return _coeffs.head<3>()                    ; }
            virtual Eigen::Matrix<Scalar,3,1> dir() const { return _coeffs.segment<3>(3); }

            //! \brief  Returns the normal, that is stored at the first three coordinates of the internal storage.
            //! \return The plane normal as a 3D Eigen::Vector Map.
            inline Eigen::Matrix<Scalar,Dim,1>::ConstFixedSegmentReturnType<3>::Type normal() const { return _coeffs.segment<3>(3); }
            //virtual Eigen::Matrix<Scalar,3,1> normal() const { return _coeffs.segment<3>(3); }

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
            static inline PlanePrimitive fromFileEntry( std::vector<Scalar> const& entries )
            {
                return PlanePrimitive( /*    pos: */ Eigen::Map<const Eigen::Matrix<Scalar,3,1> >( entries.data()  , 3 ),
                                       /* normal: */ Eigen::Map<const Eigen::Matrix<Scalar,3,1> >( entries.data()+3, 3 ) );
            }

            // ____________________GEOMETRY____________________

            //! \brief              Returns point to plane distance.
            //! \param[in] point    Point to calculate distance from.
            //! \return             Distance from point to plane.
            inline Scalar
            getDistance( Eigen::Matrix<Scalar,3,1> const& point ) const
            {
                //return this->dir().dot(point) + this->_coeffs(3);
                return (point - this->pos()).dot( this->dir() );
            }

            inline Eigen::Matrix<Scalar,3,1>
            projectPoint( Eigen::Matrix<Scalar,3,1> const& point ) const
            {
                //return point - (this->getDistance(point) * this->dir() );
                return point - (this->getDistance(point) * this->dir() );
            } // projectPoint

            /*! \brief                          Calculates the length of the plane based on the points in \p cloud, masked by \p indices and the distance from point to plane \p threshold.
             *
             *                                  The method calculates the inliers and selects the most far away point from #pos() in both directions.
             *  \tparam     _PointT             Point wrapper stored in _PointContainerT.
             *  \tparam     _PointContainerT    Type to store the points to select inliers from. Concept: std::vector<\ref GF2::PointPrimitive>. Depr: pcl::PointCloud< _PointT >::Ptr.
             *  \tparam     _IndicesContainerT  Concept: std::vector<int>.
             *  \param[out] minMax              Output container holding the four endpoints.
             *  \param[in]  cloud               Point container to look for inlier points in.
             *  \param[in]  threshold           A point in \p cloud is an inlier, if it is closer than this threshold. Usually the "scale".
             *  \param[in]  indices             Optional input to specify subset of points by indices.
             *  \return                         EXIT_SUCCESS
             */
            template <typename _PointPrimitiveT, class _IndicesContainerT, typename _PointContainerT>
            int
            getExtent( std::vector<Eigen::Matrix<Scalar,3,1> >       & minMax
                     , _PointContainerT                         const& cloud
                     , double                                   const  threshold          = 0.01
                     , _IndicesContainerT                       const* indices_arg        = NULL
                     , bool                                     const  force_axis_aligned = false ) const
            {
                typedef Eigen::Matrix<Scalar,3,1> Position;

#ifdef GF2_USE_PCL

                std::vector<int> inliers;
                {
                    inliers.reserve( cloud.size() );
                    const int stop_at = indices_arg ? indices_arg->size() : cloud.size();
                    for ( int i = 0; i != stop_at; ++i )
                    {
                        const int pid = indices_arg ? (*indices_arg)[i] : i;
                        if ( this->getDistance( cloud[pid].template pos() ) < threshold )
                            inliers.push_back( pid );
                    }

                    //std::cout << "[" << __func__ << "]: " << "indices.size(): " << stop_at << " / " << cloud.size() << std::endl;
                }

                // check size
                if ( !inliers.size() ) return EXIT_FAILURE;

                // project cloud
                _PointContainerT on_plane_cloud;
                on_plane_cloud.reserve( inliers.size() );
                for ( int pid_id = 0; pid_id != inliers.size(); ++pid_id )
                {
                    const int pid = inliers[ pid_id ];
                    on_plane_cloud.push_back( _PointPrimitiveT(this->projectPoint(cloud[pid].template pos()), cloud[pid].template dir()) );
                }
#if 0
                // debug
                {
                    pcl::PointCloud<pcl::PointXYZRGB>::Ptr c ( new pcl::PointCloud<pcl::PointXYZRGB>() );
                    pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer("on_plane_cloud") );
                    vptr->setBackgroundColor( .5, .6, .6 );
                    for ( int pid = 0; pid != on_plane_cloud.size(); ++pid )
                    {
                        pcl::PointXYZRGB pnt;
                        pnt.getVector3fMap() = on_plane_cloud[pid].template pos();
                        pnt.r = 255; pnt.g = 0; pnt.b = 0;
                        c->push_back( pnt );
                    }
                    vptr->addPointCloud( c, "onplane");
                    vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3.f, "onplane" );

                    pcl::PointCloud<pcl::PointXYZRGB>::Ptr c1 ( new pcl::PointCloud<pcl::PointXYZRGB>() );
                    for ( int pid_id = 0; pid_id != inliers.size(); ++pid_id )
                    {
                        pcl::PointXYZRGB pnt;
                        pnt.getVector3fMap() = cloud[ inliers[pid_id] ].template pos();
                        pnt.r = 0; pnt.g = 0; pnt.b = 255;
                        c1->push_back( pnt );
                    }
                    vptr->addPointCloud( c1, "origcloud" );
                    vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3.f, "origcloud" );

                    char plane_name[255];
                    sprintf( plane_name, "plane%03d", 0 );
                    vptr->addPlane( *(this->modelCoefficients()), plane_name, 0 );
                    vptr->spin();
                }
#endif

                Eigen::Matrix<Scalar,4,4> frame; // 3 major vectors as columns, and the fourth is the centroid
                {
                    processing::PCA<_IndicesContainerT>( frame, on_plane_cloud, /* indices: */ NULL ); // no indices needed, already full cloud
//                    std::cout << "frame: " << frame << std::endl;

                    if ( force_axis_aligned )
                    {
                        // get the unit axis that is most perpendicular to the 3rd dimension of the frame
                        std::pair<Position,Scalar> dim3( Position::Zero(), Scalar(FLT_MAX) );
                        {
                            Scalar tmp;
                            if ( (tmp=std::abs(frame.col(2).template head<3>().dot( Position::UnitX() ))) < dim3.second ) { dim3.first = Position::UnitX(); dim3.second = tmp; }
                            if ( (tmp=std::abs(frame.col(2).template head<3>().dot( Position::UnitY() ))) < dim3.second ) { dim3.first = Position::UnitY(); dim3.second = tmp; }
                            if ( (tmp=std::abs(frame.col(2).template head<3>().dot( Position::UnitZ() ))) < dim3.second ) { dim3.first = Position::UnitZ(); dim3.second = tmp; }
                        }
                        frame.col(0).head<3>() = frame.col(2).head<3>().cross( dim3.first             ).normalized();
                        frame.col(1).head<3>() = frame.col(2).head<3>().cross( frame.col(0).head<3>() ).normalized();
                    }
                }
//                std::cout << "frame2: " << frame << std::endl;

                _PointContainerT local_cloud;
                processing::cloud2Local<_PointPrimitiveT,_IndicesContainerT>( local_cloud, frame, on_plane_cloud, /* indices: */ NULL ); // no indices needed, it's already a full cloud
//                for ( int pid = 0; pid != local_cloud.size(); ++pid )
//                {
//                    std::cout << "on_plane_cloud[" << pid << "]: " << on_plane_cloud[pid].toString() << std::endl;
//                    std::cout << "local_cloud[" << pid << "]: " << local_cloud[pid].toString() << std::endl;
//                }

                _PointPrimitiveT min_pt, max_pt;
                processing::getMinMax3D<_IndicesContainerT>( min_pt, max_pt, local_cloud, /* indices: */ NULL );
//                std::cout << "min_pt: " << min_pt.toString()
//                          << ", max_pt: " << max_pt.toString() << std::endl;

                minMax.resize( 4 );
                minMax[0]    = minMax[1] = min_pt.template pos();
                minMax[1](1)             = max_pt.template pos()(1);
                minMax[2]    = minMax[3] = max_pt.template pos();
                minMax[3](1)             = min_pt.template pos()(1);
//                for ( int d = 0; d != minMax.size(); ++d )
//                std::cout << "minMax[" << d << "]: " << minMax[d].transpose() << std::endl;

                for ( int d = 0; d != 4; ++d )
                {
                    // to world
                    minMax[d] = (frame * (Eigen::Matrix<Scalar,4,1>() << minMax[d], Scalar(1)).finished()).template head<3>();
//                    std::cout << "minMax2[" << d << "]: " << minMax[d].transpose() << std::endl;
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

            /*! \brief Draws the Plane from extrema
             *  \tparam PointsT Concept: std::vector<pcl::PointXYZ>.
             */
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
               // v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, .9, plane_name);
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

            /*! \brief Draws plane.
             * \tparam  _PointPrimitiveT   Concept: \ref GF2::PointPrimitive.
             * \tparam _PointContainerT   Concept: std::vector< _PointPrimitiveT >.
             * \tparam _IndicesContainerT Concept: std::vector<int>.
             */
            template <class _PointPrimitiveT, class _PointContainerT, class _IndicesContainerT> static inline int
            draw( PlanePrimitive                        const& plane
                , _PointContainerT                      const& cloud
                , Scalar                                const  radius
                , _IndicesContainerT                    const* indices
                , pcl::visualization::PCLVisualizer::Ptr       v
                , std::string                           const& plane_name
                , double                                const  r
                , double                                const  g
                , double                                const  b
                , int                                   const  viewport_id = 0
                , Scalar                                const  stretch = Scalar( 1. )
                )
            {
                int err     = EXIT_SUCCESS;
                //if ( stretch != Scalar(1.) )
                //    std::cerr << "[" << __func__ << "]: " << "WARNING, Stretch for planes is unimplemented!!!" << std::endl;

                typedef Eigen::Matrix<Scalar,3,1> Position;

                std::vector<Position> minMax;
                int      it          = 0;
                int      max_it      = 10;
                Scalar  tmp_radius  = radius;
                do
                {
                    err = plane.getExtent<_PointPrimitiveT>( minMax
                                                           , cloud
                                                           , tmp_radius
                                                           , indices
                                                           , /* force_axis_aligned: */ true );
                    tmp_radius *= 2.f;
                } while ( (minMax.size() < 2) && (++it < max_it) );

                // if error or couldn't find a scale that was big enough to find "inliers"
                if ( (EXIT_SUCCESS != err) || (it >= max_it) )
                {
                    std::cerr << "[" << __func__ << "]: " << "plane.getExtent exceeded max radius increase iteration count...drawing unit " << plane.toString() << std::endl;
                    v->addPlane( *plane.modelCoefficients(), plane_name, 0 );
                }
                else
                {
                    std::vector<pcl::PointXYZ> ps;
                    for ( int i = 0; i != minMax.size(); ++i )
                    {
                        pcl::PointXYZ pnt;
                        pnt.x = minMax[i](0); pnt.y = minMax[i](1); pnt.z = minMax[i](2);
                        ps.push_back( pnt );
                    }
                    err += draw( ps, v, plane_name, r, g, b, viewport_id );
                }

                return err;
            } //...draw()
#endif //...GF2_USE_PCL

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
        _coeffs.template head<3>() = pnt;
        _coeffs.template segment<3>(3) = normal;
        //_coeffs.template segment<3>(0) = normal.normalized();
        //_coeffs                    (3) = Scalar(-1) * _coeffs.template head<3>().dot( pnt.template head<3>() ); // distance
    }

#   ifdef GF2_USE_PCL
    inline pcl::ModelCoefficients::Ptr
    PlanePrimitive::modelCoefficients() const
    {
        pcl::ModelCoefficients::Ptr model_coeffs( new pcl::ModelCoefficients() );
       model_coeffs->values.resize(Dim);
       // copy normal
       std::copy( _coeffs.data()+3, _coeffs.data()+6, model_coeffs->values.begin() );
       // calculate distance from origin
       model_coeffs->values[3] = Scalar(-1.) * this->normal().dot( this->pos() );
       if ( model_coeffs->values[3] != this->getDistance(Eigen::Matrix<Scalar,3,1>::Zero()) )
           std::cout << "these should be similar: " << this->getDistance(Eigen::Matrix<Scalar,3,1>::Zero())
                        << " , " << model_coeffs->values[3] << std::endl;
       return model_coeffs;
    }
#   endif // GF2_USE_PCL

} //...ns GF2

#endif // __GF2_INC_PLANEPRIMITIVE_HPP__

#endif // __PLANEPRIMITIVE_H__

