#ifndef RAPTER_LINEPRIMITIVE_HPP
#define RAPTER_LINEPRIMITIVE_HPP

#include "rapter/primitives/linePrimitive.h"

namespace rapter
{
    /*! \brief              Creates a LinePrimitive from a normal and a distance from origin.
     *  \tparam DerivedT    A vector of 3 scalars representing a 3D direction. Concept: Eigen::Matrix<Scalar,3,1>
     */
    template <typename DerivedT> inline int
    LinePrimitive::generateFrom( LinePrimitive                & out
                               , DerivedT                const& normal
                               , LinePrimitive::Scalar   const  distanceFromOrigin )
    {
        out.coeffs().segment<3>(3) = normal.cross( Position::UnitZ() ).normalized();
        out.coeffs().segment<3>(0) = Position::Zero() - normal * distanceFromOrigin;

        return EXIT_SUCCESS;
    } //...generatreFrom()

    /*! \brief                          Calculates the length of the line based on the points in \p cloud, masked by \p indices and the distance from point to line \p threshold.
     *
     *                                  The method calculates the inliers and selects the most far away point from #pos() in both directions.
     *  \tparam     _PointT             Type of point stored in \p cloud. Concept: \ref rapter::PointPrimitve.
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
    LinePrimitive::getExtent( LinePrimitive::ExtremaT                       & minMax
                            , _PointContainerT                         const& cloud
                            , double                                   const  threshold         /*   = 0.01*/
                            , _IndicesContainerT                       const* indices_arg       /*   = NULL*/
                            , bool                                     const  force_axis_aligned/*   = false*/ ) const
    {
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
        for ( UPidT pid_id = 0; pid_id != inliers.size(); ++pid_id )
        {
            const UPidT pid = inliers[ pid_id ];
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

    /*! \brief              Draws a line as long as it finds 2 points inside its radius. The radius is doubled 10 times until it gives up looking for more points.
     *  \tparam _PointPrimitiveT   Concept: \ref rapter::PointPrimitive.
     *  \tparam _PointContainerT   Concept: std::vector< _PointPrimitiveT >.
     *  \tparam _IndicesContainerT Concept: std::vector<int>.
     *  \param  stretch     Multiplies the final length with this number, so that the line is a bit longer than the distance between its extrema.
     */
    template <class _PointPrimitiveT, class _IndicesContainerT, class _PointContainerT>
    int
    LinePrimitive::draw( LinePrimitive                    const& line
                       , _PointContainerT                 const& cloud
                       , LinePrimitive::Scalar            const  radius
                       , _IndicesContainerT               const* indices
                       , pcl::visualization::PCLVisualizer::Ptr  v
                       , std::string                      const  plane_name
                       , double                           const  r
                       , double                           const  g
                       , double                           const  b
                       , int                              const  viewport_id
                       , LinePrimitive::Scalar            const  stretch
                       , int                              const  /*draw_mode*/
                       , float                            const  /*hull_alpha*/
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

        v->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 4.0, plane_name, 0 );

        return err;
    } //...draw()


} //...ns rapter

#endif // GO_LINEPRIMITIVE_HPP
