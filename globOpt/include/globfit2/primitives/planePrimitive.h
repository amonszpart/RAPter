#ifndef __GF2_PLANEPRIMITIVE_H__
#define __GF2_PLANEPRIMITIVE_H__

#include <Eigen/Dense>
#include "globfit2/optimization/energyFunctors.h"
#include "globfit2/primitives/primitive.h"
#include "globfit2/processing/util.hpp" // pca, getPopulationOf()

#ifdef GF2_USE_PCL
#   include "pcl/ModelCoefficients.h"
#   include "pcl/sample_consensus/sac_model_plane.h"
#   include "pcl/common/common.h"
#   include "pcl/visualization/pcl_visualizer.h"
#   include "pcl/surface/concave_hull.h"
#endif

//#include "omp.h"

namespace GF2
{
    //! \brief   Class to wrap a plane with. Implements #pos() and #dir() functions, and a constructor taking a position and a direction.
    //!
    //!          Stores 3D normal at the first three coeffs, and distance from origin at the fourth coordinate.
    //template <typename _Scalar>
    class PlanePrimitive : public ::GF2::Primitive<3,6>
    {
            typedef ::GF2::Primitive<3,6> ParentT;
            //using ParentT::_coeffs;
            //using ParentT::Scalar;
            //using ParentT::Dim;
        public:
            //enum { Dim = ParentT::Dim };
            typedef ParentT::Scalar Scalar;
            typedef std::vector<Eigen::Matrix<Scalar,3,1> > ExtentsT;

            // ____________________CONSTRUCT____________________

            inline PlanePrimitive() : ParentT() {}

            //! \brief Constructor that takes raw data in Eigen format as input.
            inline PlanePrimitive( Eigen::Matrix<Scalar,Dim,1> coeffs ) : ParentT( coeffs ) {}

            //! \brief Constructor that takes raw data in std::vector format as input.
            inline PlanePrimitive( std::vector<Scalar> const& coeffs ) : ParentT( coeffs ) {}

            /*! \brief Creates PlanePrimitive from point on plane and direction.
             *  \param[in] p0    Point on plane.
             *  \param[in] dir   Plane normal.
             */
            PlanePrimitive( Eigen::Matrix<Scalar,3,1> pnt, Eigen::Matrix<Scalar,3,1> normal );

            /*! \brief Creates plane primitive from an ordered local frame
             */
            inline PlanePrimitive( Eigen::Matrix<Scalar,3,1> const& centroid, Eigen::Matrix<Scalar,3,1> const& eigen_values, Eigen::Matrix<Scalar, 3, 3> const& eigen_vectors );

            /*! \brief Called from \ref GF2::CandidateGenerator::generate to create a new candidate from this and \p other.
             *         Create back-rotated version of the other at this position.
             *  \param[out] out              Created output primitive.
             *  \param[in]  other            Primitive, who's direction we want to use to generate something at our location.
             *  \param[in]  closest_angle_id The id of closest perfect angle between us and \p other. This we want to rotate back by.
             *  \param[in]  angles           List of desired angles, one of its entries is referenced by closest_angle_id.
             *  \return \p out is only valid, if return true. We might decide, that there's nothing to create, in which case we return false.
             */
            template <class _AngleContainerT>
            bool generateFrom( PlanePrimitive         & out
                             , PlanePrimitive    const& other
                             , int               const  closest_angle_id
                             , _AngleContainerT  const& angles
                             , Scalar            const  angle_multiplier = Scalar(1.)
                             ) const;
            template <typename DerivedT> static
            int generateFrom( PlanePrimitive                 & out
                            , DerivedT                  const& normal
                            , Scalar                    const  distanceFromOrigin );

            // ____________________VIRTUALS____________________
            /*! \brief  Compulsory virtual overload of position getter. The position of the plane is calculated on the fly from the formula N . x0 + d = 0.
             *  \return The position of the plane as a 3D Eigen::Vector.
             */
            //virtual Eigen::Matrix<Scalar,3,1> pos() const { return _coeffs.template head<3>(); }
            inline typename Eigen::Matrix<Scalar,Dim,1>::ConstFixedSegmentReturnType<3>::Type pos() const { return _coeffs.template head<3>(); }

            /*! \brief  Compulsory virtual overload of orientation getter. The orientation of the plane is the normal stored at the first three coordinates of #_coeffs.
             *  \return The normal of the plane as a 3D Eigen::Vector.
             */
            //virtual Eigen::Matrix<Scalar,3,1> dir() const { return _coeffs.template segment<3>(3); }
            inline typename Eigen::Matrix<Scalar,Dim,1>::ConstFixedSegmentReturnType<3>::Type dir() const { return _coeffs.template segment<3>(3); }

            /*! \brief  Returns the normal, that is stored at the first three coordinates of the internal storage.
             *  \return The plane normal as a 3D Eigen::Vector Map.
             */
            inline typename Eigen::Matrix<Scalar,Dim,1>::ConstFixedSegmentReturnType<3>::Type normal() const { return _coeffs.segment<3>(3); }

            // _______________________IO_______________________
            /*! \brief Used in \ref io::readPrimitives to determine how many floats to parse from one entry.
             */
            static inline int getFileEntryLength() { return 6; }

            /*! \brief Output <x0,n>, location and normal for the plane to a string that does *not* have an endline at the end.
             */
            std::string toFileEntry() const;

            /*! \brief              Used in \ref io::readPrimitives to create the primitive from the read floats.
             *  \param[in] entries  Contains <x0,n> to create the plane from.
             */
            static PlanePrimitive fromFileEntry( std::vector<Scalar> const& entries );

            // ____________________GEOMETRY____________________

            /*! \brief              Returns point to plane distance.
             *  \param[in] point    Point to calculate distance from.
             *  \return             Distance from point to plane.
             */
            inline Scalar
            getDistance( Eigen::Matrix<Scalar,3,1> const& point ) const {  return (point - this->pos()).dot( this->dir() ); }

            /*! \brief              Helps calculating plane to plane distance.
             *  \param[in] extrema  Extrema of this primitive.
             *  \param[in] pnt      One of the extrema of the other primitive.
             *  \return             Distance from point to plane.
             */
            Scalar
            getFiniteDistance( ExtentsT const& extrema, Position const& pnt ) const;

            int to4Coeffs( std::vector<Scalar> &coeffs ) const;

            Eigen::Matrix<Scalar,3,1>
            projectPoint( Eigen::Matrix<Scalar,3,1> const& point ) const;

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
            getExtent( ExtentsT                                      & minMax
                     , _PointContainerT                         const& cloud
                     , double                                   const  threshold          = 0.01
                     , _IndicesContainerT                       const* indices_arg        = NULL
                     , bool                                     const  force_axis_aligned = false ) const;

            /*! \brief Calculates size, a bit smarter, than taking the area of #getExtent().
             *  \tparam MatrixDerived   Concept: Eigen::Matrix<_Scalar,-1,1>.
             */
            template <class _IndicesContainerT, typename MatrixDerived, typename _Scalar, class _PointContainerT >
            MatrixDerived& getSpatialSignificance( MatrixDerived             & in
                                                 , _PointContainerT     const& points
                                                 , _Scalar             const /*scale*/
                                                 , _IndicesContainerT        * indices        = NULL
                                                 , bool                 const  return_squared = false ) const;

            // ____________________DRAWING____________________
#ifdef GF2_USE_PCL
            inline pcl::ModelCoefficients::Ptr
            modelCoefficients() const;

            /*! \brief Draws the Plane from extrema
             *  \tparam PointsT Concept: std::vector<pcl::PointXYZ>.
             */
            template <class PointsT> static inline int
            draw( PointsT                                   const& corners
                  , pcl::visualization::PCLVisualizer::Ptr         v
                  , std::string                             const  plane_name
                  , double                                  const  r
                  , double                                  const  g
                  , double                                  const  b
                  , int                                     const  viewport_id = 0
                  );

            /*! \brief Extract convec hull and display
             *  \tparam PclCloudT          Concept: pcl::PointCloud<pcl::PointXYZ>::Ptr plane_polygon_cloud_ptr( new pcl::PointCloud<pcl::PointXYZ> ).
             *  \tparam _PointContainerT   Concept: std::vector<pcl::PointXYZ>.
             *  \tparam _IndicesContainerT Concept: std::vector<int>.
             */
            template < class PclCloudT, class _PointContainerT, class _IndicesContainerT> static inline int
            getHull( PclCloudT                & plane_polygon_cloud
                   , PlanePrimitive      const& plane
                   , _PointContainerT    const& points
                   , _IndicesContainerT  const* indices
                   , float               const  alpha = 2.f
                   , pcl::PolygonMesh         * out_mesh = NULL
                   );

            /*! \brief Extract convec hull and display
             *  \tparam PointsT Concept: std::vector<pcl::PointXYZ>.
             */
            template <class _PointContainerT, class _IndicesContainerT> static inline int
            drawConvex( pcl::visualization::PCLVisualizer::Ptr v
                      , PlanePrimitive      const& plane
                      , _PointContainerT    const& cloud
                      , _IndicesContainerT  const* indices
                      , std::string         const  plane_name
                      , double              const  r
                      , double              const  g
                      , double              const  b
                      , int                 const  viewport_id = 0
                      , float               const  alpha       = 2.f
                      );

            /*! \brief Draws plane.
             * \tparam  _PointPrimitiveT   Concept: \ref GF2::PointPrimitive.
             * \tparam _PointContainerT   Concept: std::vector< _PointPrimitiveT >.
             * \tparam _IndicesContainerT Concept: std::vector<int>.
             */
            template <class _PointPrimitiveT, class _PointContainerT, class _IndicesContainerT> static int
            draw( PlanePrimitive                        const& plane
                , _PointContainerT                      const& cloud
                , Scalar                                const  radius
                , _IndicesContainerT                    const* indices
                , pcl::visualization::PCLVisualizer::Ptr       v
                , std::string                           const& plane_name
                , double                                const  r
                , double                                const  g
                , double                                const  b
                , int                                   const  viewport_id  = 0
                , Scalar                                const  stretch      = Scalar( 1. )
                , int                                   const  draw_mode    = 0
                , Scalar                                const  alpha        = 2.
                );
#endif //...GF2_USE_PCL

        bool gidUnset() const;
    }; //...class PlanePrimitive
} //...ns GF2



// ________________________________________________________HPP_________________________________________________________

#if 0
#ifndef __GF2_INC_PLANEPRIMITIVE_HPP__
#define __GF2_INC_PLANEPRIMITIVE_HPP__
namespace GF2
{
    PlanePrimitive::PlanePrimitive( Eigen::Matrix<Scalar, 3, 1> pnt, Eigen::Matrix<Scalar, 3, 1> normal )
    {
        _coeffs.template head<3>() = pnt;
        _coeffs.template segment<3>(3) = normal.normalized();
        //_coeffs.template segment<3>(0) = normal.normalized();
        //_coeffs                    (3) = Scalar(-1) * _coeffs.template head<3>().dot( pnt.template head<3>() ); // distance
    }

    /*! \brief Creates plane primitive from an ordered local frame
     */
    PlanePrimitive::PlanePrimitive( Eigen::Matrix<Scalar,3,1> const& centroid, Eigen::Matrix<Scalar,3,1> const& eigen_values, Eigen::Matrix<Scalar, 3, 3> const& eigen_vectors )
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
    template <class _AngleContainerT> bool
    PlanePrimitive::generateFrom( PlanePrimitive         & out
                                , PlanePrimitive    const& other
                                , int               const  closest_angle_id
                                , _AngleContainerT  const& angles
                                , Scalar            const  angle_multiplier /* = Scalar(1.) */
                                ) const
    {
        // if not 0 or M_PI, meaning not parallel
        if ( (closest_angle_id != 0) && (closest_angle_id != angles.size()-1) )
        {
            //std::cerr << "[" << __func__ << "]: " << "rotating plane by angle " << angles[closest_angle_id] << std::endl;
            //return false;
            //Scalar const angle = angles[ closest_angle_id ];
            Scalar angle = angles[ closest_angle_id ];
#warning This has to be tested for angle sign
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
        out.copyTagsFrom( *this );
        // copy position id from self
        //out.setTag( GID, this->getTag(GID) );
        // copy direction id from the other
        out.setTag( TAGS::DIR_GID, other.getTag(TAGS::DIR_GID) );
        // erase chosen tag - this is a new candidate
        out.setTag( TAGS::STATUS, STATUS_VALUES::UNSET );

        return true;
    } //...generateFrom

    int PlanePrimitive::generateFrom( PlanePrimitive                 & out
                                    , Eigen::Matrix<Scalar,3,1> const& normal
                                    , Scalar                    const  distanceFromOrigin )
    {
        out.coeffs().segment<3>(3) = normal.normalized();
        out.coeffs().segment<3>(0) = Eigen::Matrix<Scalar,3,1>::Zero() - normal * distanceFromOrigin;
        return EXIT_SUCCESS;
    } //...generatreFrom()

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
    PlanePrimitive::getExtent( ExtentsT                                      & minMax
                             , _PointContainerT                         const& cloud
                             , double                                   const  threshold          /*= 0.01*/
                             , _IndicesContainerT                       const* indices_arg        /*= NULL*/
                             , bool                                     const  force_axis_aligned /*= false */) const
    {
        typedef Eigen::Matrix<Scalar,3,1> Position;

        if ( this->_extents.isUpdated() )
        {
            minMax = _extents.get();
            return 0;
        }

#ifdef GF2_USE_PCL

        std::vector<LidT> inliers;
        {
            inliers.reserve( cloud.size() );
            const PidT stop_at = indices_arg ? indices_arg->size() : cloud.size();
            for ( PidT i = 0; i != stop_at; ++i )
            {
                const PidT pid = indices_arg ? (*indices_arg)[i] : i;
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
        for ( PidT pid_id = 0; pid_id != inliers.size(); ++pid_id )
        {
            const PidT pid = inliers[ pid_id ];
            on_plane_cloud.push_back( _PointPrimitiveT(this->projectPoint(cloud[pid].template pos()), cloud[pid].template dir()) );
        }

#if 0
        // debug
        {
            pcl::PointCloud<pcl::PointXYZRGB>::Ptr c ( new pcl::PointCloud<pcl::PointXYZRGB>() );
            pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer("on_plane_cloud") );
            vptr->setBackgroundColor( .5, .6, .6 );
            for ( PidT pid = 0; pid != on_plane_cloud.size(); ++pid )
            {
                pcl::PointXYZRGB pnt;
                pnt.getVector3fMap() = on_plane_cloud[pid].template pos();
                pnt.r = 255; pnt.g = 0; pnt.b = 0;
                c->push_back( pnt );
            }
            vptr->addPointCloud( c, "onplane");
            vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3.f, "onplane" );

            pcl::PointCloud<pcl::PointXYZRGB>::Ptr c1 ( new pcl::PointCloud<pcl::PointXYZRGB>() );
            for ( PidT pid_id = 0; pid_id != inliers.size(); ++pid_id )
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

#if 0
        Scalar step = Scalar(1. * M_PI) / Scalar(180.);

        Scalar limits[2] = { Scalar(0.), M_PI/Scalar(2.) };
        for ( int it = 0; it != 2; ++it, step /= Scalar(5.) )
        {
            std::pair <Scalar,Scalar> min_volume; // <ang, volume>
            min_volume.first  = Scalar(-1.);
            min_volume.second = Scalar(FLT_MAX);
            int i = 0;
            for ( Scalar ang = limits[0]; ang < limits[1]; ang += step, ++i )
            {
                std::cout << "rot by " << ang << ", step: " << step << std::endl;
                // rotated frame
                Eigen::Matrix4f tmp_frame = frame;

                // rotate frame around up_vector by ang
                Eigen::AngleAxisf rot( ang, tmp_frame.block<3,1>(0,2) );
                tmp_frame.block<3,1>(0,0) = rot * tmp_frame.block<3,1>(0,0);
                tmp_frame.block<3,1>(0,1) = rot * tmp_frame.block<3,1>(0,1);

                // calculate volume
                _PointContainerT local_cloud;
                processing::cloud2Local<_PointPrimitiveT,_IndicesContainerT>( local_cloud, frame, on_plane_cloud, /* indices: */ NULL ); // no indices needed, it's already a full cloud

                _PointPrimitiveT min_pt, max_pt;
                processing::getMinMax3D<_IndicesContainerT>( min_pt, max_pt, local_cloud, /* indices: */ NULL );

                Eigen::Vector3f diag    = max_pt.template pos() - min_pt.template pos();

                //get location of minimum
                Eigen::MatrixXf::Index minRow, minCol;
                diag.minCoeff( &minRow, &minCol );
                std::cout << "minRow: " << minRow << std::endl;
                float volume = 0.;
                if ( minRow == 1 )
                    volume = diag(0) * diag(2);
                else if ( minRow )
                    volume = diag(0) * diag(1);
                else
                    volume = diag(1) * diag(2);
                // select min
                if ( volume < min_volume.second )
                {
                    min_volume.first  = ang;
                    min_volume.second = volume;
                }
            }

//            if ( true )
//                std::cout << "min_volume.first (angle): " << min_volume.first << " radians "
//                          << min_volume.first * Scalar(180.)/Scalar(M_PI) << " degrees, "  << min_volume.second << " volume" << std::endl;

            // selected apply rotation
            Eigen::AngleAxisf rot( min_volume.first, frame.block<3,1>(0,2) );
            frame.block<3,1>(0,0) = rot * frame.block<3,1>(0,0);
            frame.block<3,1>(0,1) = rot * frame.block<3,1>(0,1);

            // modify lookup around chosen angle for next iteration
            limits[0] = -step/2.f;
            limits[1] = limits[0] + step;
        }
#endif

        _PointContainerT local_cloud;
        processing::cloud2Local<_PointPrimitiveT,_IndicesContainerT>( local_cloud, frame, on_plane_cloud, /* indices: */ NULL ); // no indices needed, it's already a full cloud

        _PointPrimitiveT min_pt, max_pt;
        processing::getMinMax3D<_IndicesContainerT>( min_pt, max_pt, local_cloud, /* indices: */ NULL );

        minMax.resize( 4 );
        minMax[0]    = minMax[1] = min_pt.template pos();
        minMax[1](1)             = max_pt.template pos()(1);
        minMax[2]    = minMax[3] = max_pt.template pos();
        minMax[3](1)             = min_pt.template pos()(1);

        for ( int d = 0; d != 4; ++d )
        {
            // to world
            minMax[d] = (frame * (Eigen::Matrix<Scalar,4,1>() << minMax[d], Scalar(1)).finished()).template head<3>();
        }

        this->_extents.update( minMax );

        return EXIT_SUCCESS;
#           else
        std::cerr << "[" << __func__ << "]: " << "Needs PCL to work!!!!"
        return exit_FAILURE;
#           endif //...GF2_USE_PCL
    } //...getExtent()

    /*! \brief Calculates size, a bit smarter, than taking the area of #getExtent().
     *  \tparam MatrixDerived   Concept: Eigen::Matrix<_Scalar,-1,1>.
     */
    template <class _IndicesContainerT, typename MatrixDerived, typename _Scalar, class _PointContainerT > MatrixDerived&
    PlanePrimitive::getSpatialSignificance( MatrixDerived            & in
                                          , _PointContainerT    const& points
                                          , _Scalar             const /*scale*/
                                          , _IndicesContainerT       * indices       /*= NULL*/
                                          , bool                const return_squared /*= false */) const
    {
#warning "[PlanePrimitive::getSpatialSignificance] TODO: Project onto normals."
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

    /*! \brief Draws the Plane from extrema
     *  \tparam PointsT Concept: std::vector<pcl::PointXYZ>.
     */
    template <class PointsT> int
    PlanePrimitive::draw( PointsT                               const&  corners
                        , pcl::visualization::PCLVisualizer::Ptr        v
                        , std::string                           const   plane_name
                        , double                                const   r
                        , double                                const   g
                        , double                                const   b
                        , int                                   const   viewport_id /* = 0 */
                        )
    {
#   ifdef GF2_USE_PCL
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

//        #pragma omp critical PCLVIS
        {
            // draw plane polygon
            v->addPolygon<pcl::PointXYZ>( plane_polygon_cloud_ptr, r,g,b, plane_name, viewport_id );
            // v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, .9, plane_name);
            v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_SURFACE, plane_name );
        }

        // show normal
        /* if ( drawNormal )
        {
            v->addArrow( am::util::pcl::asPointXYZ( centroid + plane.Normal() * .2f ), // end point
                         am::util::pcl::asPointXYZ( centroid                        ), // start point
                         r*3., g*3., b*3.,
                         false,                                                                   // show length
                         std::string(title) + "_normal"                                                  // cloud id
                         , viewport_id );
        }

        // draw coorindate system
        if ( drawFrame )
        {
            for ( int axis_id = 0; axis_id < 3; ++axis_id )
            {
                char arrow_name[1024];
                sprintf( arrow_name, "_arrow%d", axis_id );

                v->addArrow( am::util::pcl::asPointXYZ( plane.Frame().block<3,1>(0,3) + plane.Frame().block<3,1>(0,axis_id).normalized() *.2f ), // end point
                             am::util::pcl::asPointXYZ( plane.Frame().block<3,1>(0,3)                             ), // start point
                             (2-axis_id  )%3 *.45 + .1,                                               // r
                             (2-axis_id+1)%3 *.45 + .1,                                               // g
                             (2-axis_id+2)%3 *.45 + .1,                                               // b
                             false,                                                                   // show length
                             std::string(title) + arrow_name                                          // cloud id
                             , viewport_id );
                //v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_, 3., std::string(title) + arrow_name );
            }
        } */

        return EXIT_SUCCESS;
#   else
        return EXIT_FAILURE;
#   endif // GF2_USE_PCL
    } //...draw()

    /*! \brief Draws plane.
     * \tparam  _PointPrimitiveT   Concept: \ref GF2::PointPrimitive.
     * \tparam _PointContainerT   Concept: std::vector< _PointPrimitiveT >.
     * \tparam _IndicesContainerT Concept: std::vector<int>.
     */
    template <class _PointPrimitiveT, class _PointContainerT, class _IndicesContainerT> int
    PlanePrimitive::draw( PlanePrimitive                        const& plane
                        , _PointContainerT                      const& cloud
                        , Scalar                                const  radius
                        , _IndicesContainerT                    const* indices
                        , pcl::visualization::PCLVisualizer::Ptr       v
                        , std::string                           const& plane_name
                        , double                                const  r
                        , double                                const  g
                        , double                                const  b
                        , int                                   const  viewport_id  /* = 0*/
                        , Scalar                                const  stretch      /* = Scalar( 1. ) */
                        , int                                   const  draw_mode    /* = 0*/
                        , Scalar                                const  alpha        /* = 2.*/
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
                                                   , /* force_axis_aligned: */ (draw_mode == 1) /*true*/ );
            tmp_radius *= 2.f;
        } while ( (minMax.size() < 2) && (++it < max_it) );

        // if error or couldn't find a scale that was big enough to find "inliers"
        if ( (EXIT_SUCCESS != err) || (it >= max_it) )
        {
            std::cerr << "[" << __func__ << "]: " << "plane.getExtent exceeded max radius increase iteration count...drawing unit " << plane.toString() << std::endl;
            v->addPlane( *plane.modelCoefficients(), plane.pos()(0), plane.pos()(1), plane.pos()(2), plane_name, 0 );
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
            if ( draw_mode <= 1 ) // 0: classic, 1: classic axis_aligned, 2: qhull
                err += draw( ps, v, plane_name, r, g, b, viewport_id );
            else if ( draw_mode == 2 )
                err += drawConvex( v, plane, cloud, indices, plane_name, r, g, b, viewport_id, alpha );
            else
            {
                std::cerr << "can't recognize draw mode " << draw_mode << std::endl;
                err = EXIT_FAILURE;
            }
        }

        //v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_SHADING, pcl::visualization::PCL_VISUALIZER_SHADING_FLAT, plane_name );

        return err;
    } //...draw()

    /*! \brief Extract convec hull and display
     *  \tparam PclCloudT Concept: pcl::PointCloud<pcl::PointXYZ>::Ptr plane_polygon_cloud_ptr( new pcl::PointCloud<pcl::PointXYZ> ).
     *  \tparam PointsT Concept: std::vector<pcl::PointXYZ>.
     */
    template < class PclCloudT, class _PointContainerT, class _IndicesContainerT> int
    PlanePrimitive::getHull( PclCloudT                & plane_polygon_cloud
                           , PlanePrimitive      const& plane
                           , _PointContainerT    const& points
                           , _IndicesContainerT  const* indices
                           , float               const  alpha    /* = 2.f  */
                           , pcl::PolygonMesh         * out_mesh /* = NULL */
                           )
    {
        typedef typename PclCloudT::PointType PclPointT;

        pcl::ConcaveHull<PclPointT>              concave_hull;                                   // object
        //typename pcl::PointCloud<PclPointT>::Ptr cloud_hull( new pcl::PointCloud<PclPointT>() );
        typename pcl::PointCloud<PclPointT> cloud_hull;
        typename pcl::PointCloud<PclPointT>::Ptr cloud_projected( new typename pcl::PointCloud<PclPointT>() );
        std::vector<pcl::Vertices>          polygons;                               // output list indexing the points from cloud_hull, in 2D this is size 1
        //pcl::PointCloud<PclPointT>::Ptr plane_polygon_cloud_ptr( new pcl::PointCloud<PclPointT> );

        cloud_projected->resize(indices->size());

        // get assigned points, project them to the plane and store as PCL cloud
        PidT i = 0;
        for ( typename _IndicesContainerT::const_iterator it = indices->begin(); it != indices->end(); ++i, ++it )
        {
            cloud_projected->at(i).getVector3fMap() = plane.projectPoint(points[*it].pos()).template cast<float>();
        }

        concave_hull.setAlpha( alpha );
        concave_hull.setInputCloud( cloud_projected );
        concave_hull.reconstruct( cloud_hull, polygons );

        if ( polygons.size() )
        {
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
        }
        else
            std::cout << "no hull" << std::endl;

        return plane_polygon_cloud.size();
    }

    /*! \brief Extract convec hull and display
     *  \tparam PointsT Concept: std::vector<pcl::PointXYZ>.
     */
    template <class _PointContainerT, class _IndicesContainerT> int
    PlanePrimitive::drawConvex( pcl::visualization::PCLVisualizer::Ptr v
                              , PlanePrimitive      const& plane
                              , _PointContainerT    const& cloud
                              , _IndicesContainerT  const* indices
                              , std::string         const  plane_name
                              , double              const  r
                              , double              const  g
                              , double              const  b
                              , int                 const  viewport_id /* = 0   */
                              , float               const  alpha       /* = 2.f */
                              )
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr plane_polygon_cloud_ptr( new pcl::PointCloud<pcl::PointXYZ> );
        if ( getHull(*plane_polygon_cloud_ptr, plane, cloud, indices, alpha) )
        {
//            #pragma omp critical PCLVIS
            {
                // draw plane polygon
                v->addPolygon<pcl::PointXYZ>( plane_polygon_cloud_ptr, r,g,b, plane_name, viewport_id );

                // v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, .9, plane_name);
                v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_SURFACE, plane_name );
            }
        }
        else
        {
            std::cerr << "[" << __func__ << "]: " << "qhull returned 0 points for plane " << plane.toString() << "..." << std::endl;
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    }

    inline int PlanePrimitive::to4Coeffs( std::vector<Scalar> &coeffs ) const
    {
        coeffs.resize( 4, Scalar(0.) );
        std::copy( _coeffs.data()+3, _coeffs.data()+6, coeffs.begin() );
        // calculate distance from origin
        coeffs[3] = Scalar(-1.) * this->normal().dot( this->pos() );
        return EXIT_SUCCESS;
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

    inline bool PlanePrimitive::gidUnset() const { return this->getTag(PlanePrimitive::TAGS::GID) == LONG_VALUES::UNSET; }
} //...ns GF2

#endif // __GF2_INC_PLANEPRIMITIVE_HPP__

#endif

#endif // __PLANEPRIMITIVE_H__

