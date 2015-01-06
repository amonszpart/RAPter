#ifndef __GF2_VISUALIZER_HPP__
#define __GF2_VISUALIZER_HPP__

#include "globfit2/globOpt_types.h"

#if GF2_USE_PCL
#   include "pcl/point_types.h"
#   include "pcl/visualization/pcl_visualizer.h"
#endif

#include "globfit2/visualization/visualization.h" // MyVisPtr
#include "globfit2/processing/util.hpp"           // getPopulations()

#include "qcqpcpp/optProblem.h"

namespace GF2 {
    //! \brief Visualizer class to show points, primitives and their relations
    //! \tparam PrimitiveContainerT     Patch-wise grouped primitive container type. Concept: vector< vector< GF2::LinePrimitive2> >.
    //! \tparam PointContainerT         Cloud of possibly oriented 3D points. Concept: vector< GF2::PointPrimitive >.
    template <class PrimitiveContainerT, class PointContainerT>
    class Visualizer
    {
            typedef pcl::PointXYZRGB            MyPCLPoint;
            typedef pcl::PointCloud<MyPCLPoint> MyPCLCloud;

        public:
            typedef Eigen::Vector3f             Colour;

            // tested draw_modes: 0, 1, 2, 4 (REPROJECT), 12 (REPROJECT | HIDE_PRIMITIVES), 28 ( REPROJECT | HIDE_PRIMITIVES | HIDE_UNASSIGNED_PTS )
            enum DRAW_MODE { SIMPLE = 0, AXIS_ALIGNED = 1, QHULL = 2, REPROJECT = 4
                           , HIDE_PRIMITIVES = 8, HIDE_UNASSIGNED_PTS = 16 };
            /*! \brief                          Visualize lines and points
             *  \tparam _Scalar                 Floating point type used in points to store data. Concept: typename PointContainerT::value_type::Scalar.
             *  \param spin                     Halts execution and allows interactive viewing, if true
             *  \param draw_cons                Draw perfect angles
             *  \param show_ids                 Draw line ids (lid,lid1)
             *  \param use_tags                 Restrict line extent to GID tagged points. case 1: colour coded points with ellipses. case 2: colour coded points with no ellipses
             *  \param[in] perfect_angle_limit  When to display gray line for "perfect angle".
             *  \param[in] print_perf_angles    Show the angles in degrees, if below perfect_angle_limit.
             *  \param[in] stretch              Elong primitives beyond their extrema by multiplying their dimensions by this number (1 == don't elong, 1.2 == elong a bit)
             *  \param[in] draw_mode            Mode0: classic, Mode1: classic, axis aligned, Mode2: qhull
             *  \param[in] hull_alpha           Alpha parameter for convex hull calculation
             *  \return             The visualizer for further display and manipulation
             */
            template <typename _Scalar> static inline vis::MyVisPtr
            show( PrimitiveContainerT  const& primitives
                , PointContainerT      const& points
                , _Scalar              const  scale
                , Colour               const& colour                = (Colour() << 0.f, 0.f, 1.f).finished()
                , Colour               const& bg_colour             = (Colour() << .1f, .1f, .1f).finished()
                , bool                 const  spin                  = true
                , std::vector<_Scalar> const* angles                = NULL
                , bool                 const  show_ids              = false
                , char                 const  use_tags              = false
                , int                  const  pop_limit             = 10
                , std::string          const& title                 = ""
                , bool                 const  show_pids             = false
                , int                  const  show_normals          = 0
                , bool                 const  show_pop              = false
                , _Scalar              const  perfect_angle_limit   = 10.e-6
                , bool                 const  print_perf_angles     = false
                , bool                 const  dir_colours           = false
                , bool                 const  no_points             = false
                , std::set<int>        const* filter_gids           = NULL
                , std::set<int>        const* filter_dids           = NULL
                , std::set<int>        const* filter_status         = NULL
                , _Scalar              const  stretch               = 1.
                , int                  const  draw_mode             = DRAW_MODE::SIMPLE
                , bool                 const  no_scale_sphere       = false
                , _Scalar              const  hull_alpha            = 2.
                , bool                 const  save_poly             = false
                , _Scalar              const  point_size            = 6.0
                , bool                 const  skip_empty            = false
                , bool                 const  show_spatial          = false
                , std::string          const  problemPath           = ""
                );

            //! \brief Shows a polygon that approximates the bounding ellipse of a cluster
            template <typename _Scalar> static inline int
            drawEllipse( vis::MyVisPtr                             vptr
                       , MyPCLCloud::Ptr                              cloud
                       , std::vector<PidT>                  const& indices
                       , _Scalar                            const  scale
                       , int                                const  prim_tag
                       , Eigen::Matrix<_Scalar,3,1>         const& prim_colour
                       );
    }; //...class Visualizer
} // ns GF2

//_________________________________________________________________
//______________________________HPP________________________________
//_________________________________________________________________

#if GF2_USE_PCL
#   include <pcl/common/common.h>        // getMinMax3D
#   include <pcl/visualization/pcl_visualizer.h>
#   include <pcl/common/pca.h>
#   include "globfit2/util/pcl_util.hpp" // pclutil::asPointXYZ
#   include "pcl/PolygonMesh.h"
#endif

#include "globfit2/my_types.h"
#include "globfit2/util/util.hpp" // util::nColoursEigen
#include "globfit2/optimization/energyFunctors.h"

namespace GF2
{
    template <class PrimitiveContainerT, class PointContainerT>
    template <typename _Scalar>
    vis::MyVisPtr
    Visualizer<PrimitiveContainerT,PointContainerT>::show( PrimitiveContainerT    const& primitives
                                                           , PointContainerT      const& points
                                                           , _Scalar              const  scale
                                                           , Colour               const& colour              /* = {0,0,1} */
                                                           , Colour               const& bg_colour
                                                           , bool                 const  spin                /* = true */
                                                           , std::vector<_Scalar> const* angles              /* = NULL */
                                                           , bool                 const  show_ids            /* = false */
                                                           , char                 const  use_tags            /* = false */
                                                           , int                  const  pop_limit           /* = 10 */
                                                           , std::string          const& title               /* = "" */
                                                           , bool                 const  show_pids           /* = false */
                                                           , int                  const  show_normals        /* = 0 */
                                                           , bool                 const  show_pop            /* = false */
                                                           , _Scalar              const  angle_limit /* = 0872664625 = 5deg*/
                                                           , bool                 const  print_perf_angles   /* = false */
                                                           , bool                 const  dir_colours         /* = false */
                                                           , bool                 const  hide_points         /* = false */
                                                           , std::set<int>        const* filter_gids         /* = NULL */
                                                           , std::set<int>        const* filter_dids         /* = NULL */
                                                           , std::set<int>        const* filter_status       /* = NULL */
                                                           , _Scalar              const  stretch             /* = 1. */
                                                           , int                  const  draw_mode           /* = 0 */
                                                           , bool                 const  no_scale_sphere     /* = false */
                                                           , _Scalar              const  hull_alpha          /* = 2. */
                                                           , bool                 const  save_poly           /* = false */
                                                           , _Scalar              const  point_size          /* = 6.0 */
                                                           , bool                 const  skip_empty          /* = false */
                                                           , bool                 const  show_spatial        /* = false */
                                                           , std::string          const  problemPath
                                                           )
    {
        bool save_hough = true;

        // TYPEDEFS
        typedef typename PrimitiveContainerT::value_type::value_type    PrimitiveT;
        typedef typename PointContainerT::value_type                    PointPrimitiveT;
        typedef          Eigen::Matrix<_Scalar,3,1>                     Position;
        typedef          std::pair< int, int     >                      LidLid1;
        typedef          std::map < int, LidLid1 >                      GidLidLid1Map;

        // --------------------------------------------------------------------
        // CONSTs
        const Eigen::Matrix<double,3,1> gray          ( (Eigen::Matrix<double,3,1>()<<.6,.6,.5).finished() );
        const _Scalar                   deg_multiplier( _Scalar(180.) / M_PI );
        const _Scalar                   text_size = scale/0.05 * 0.015;

        std::cout << "[" << __func__ << "]: " << "scale: " << scale
                  << std::endl;

        // --------------------------------------------------------------------

        // every direction or every patch gets its own colour
        const int primColourTag = dir_colours ? PrimitiveT::TAGS::DIR_GID : PrimitiveT::TAGS::GID;
        // primitives are prepared for draw_modes: SIMPLE, AXIS_ALIGNED, QHULL, so we need to remove the rest of the flags
        const int old_draw_mode = draw_mode & (SIMPLE | AXIS_ALIGNED | QHULL);
        std::cout << "old_draw_mode: " << old_draw_mode << std::endl;

        // --------------------------------------------------------------------

        GidT max_gid     = 0;
        DidT max_dir_gid = 0;
        LidT nPrimitives = 0, nbColour = 0;
        GidLidLid1Map gid2lidLid1;    // maps gid to a primitve (assuming only one per patch...)
        std::map<int, int> id2ColId;       // maps gid to a colour id

        // get the highest gid from points OR primitives, store into max_gid
        {
            for ( size_t pid = 0; pid != points.size(); ++pid )
                max_gid = std::max( max_gid, points[pid].getTag(PointPrimitiveT::TAGS::GID) );

            for ( size_t lid = 0; lid != primitives.size(); ++lid )
            {
                for ( size_t lid1 = 0; lid1 != primitives[lid].size(); ++lid1 )
                {
                    if ( primitives[lid][lid1].getTag(PrimitiveT::TAGS::STATUS) != PrimitiveT::STATUS_VALUES::SMALL ) // TODO: use filter status
                    {
                        max_gid     = std::max( max_gid    , primitives[lid][lid1].getTag(PrimitiveT::TAGS::GID    ) );
                        max_dir_gid = std::max( max_dir_gid, primitives[lid][lid1].getTag(PrimitiveT::TAGS::DIR_GID) );

                        // this id has not been referenced to a color previously
                        if ( id2ColId.find(primitives[lid][lid1].getTag(primColourTag)) == id2ColId.end())
                        {
                            id2ColId[ primitives[lid][lid1].getTag(primColourTag) ] = nbColour;
                            ++nbColour;
                        }
                        gid2lidLid1[ primitives[lid][lid1].getTag(PrimitiveT::TAGS::GID) ] = std::pair<int,int>(lid,lid1);

                        ++nPrimitives;
                    }
                }
            }
        } // max_gid, max_dir_gid, id2ColId, gid2lidLid1

        // --------------------------------------------------------------------

        // log
        {
            std::cout << "[" << __func__ << "]: "
                      << "points: "         << points.size()
                      << ", primitives: "   << nPrimitives
                      << ", max_gid: "      << max_gid
                      << ", max_dir_gid: "  << max_dir_gid
                      << ", pop-limit: "    << pop_limit
                      << ", filter-gids.size(): " << (filter_gids ? filter_gids->size() : 0)
                      << std::endl;
        }

        // --------------------------------------------------------------------

        // Check if we can use color palette
        const int paletteRequiredSize = id2ColId.size() + 1; // set this to 0, if you only want 7 colours
        // Defaults to true, since util replicates colours
        bool usePalette = id2ColId.size() <= 9; //util::paletteLightColoursEigen(paletteRequiredSize).size();
        std::vector<Eigen::Vector3f> pointColours, primColours;
        Eigen::Vector3f              unusedPointColour, unusedPrimColour;
        // Choose, and fill colour palette
        {
            if ( usePalette )
            {
                cout << "Switching to nicer color palette" << endl;

                pointColours = util::paletteMediumColoursEigen(paletteRequiredSize);
                primColours  = util::paletteDarkColoursEigen(paletteRequiredSize);

                unusedPointColour = util::paletteLightNeutralColour();
                unusedPrimColour  = util::paletteDarkNeutralColour();

            }
            else
            {
                cout << " NOT Switching to nicer color palette, need " << id2ColId.size() << " colours" << endl;
                pointColours = util::nColoursEigen( /*           count: */ nbColour
                                                    , /*         scale: */ 255.f
                                                    , /*       shuffle: */ true
                                                    , /* min_hsv_value: */ 70.f );
                primColours = pointColours;
                for ( size_t cid = 0; cid != pointColours.size(); ++cid )
                {
                    pointColours[cid](0) = std::min( pointColours[cid](0) * 1.6, 255.);
                    pointColours[cid](1) = std::min( pointColours[cid](1) * 1.6, 255.);
                    pointColours[cid](2) = std::min( pointColours[cid](2) * 1.6, 255.);
                }

                unusedPointColour = util::paletteLightNeutralColour();
                unusedPrimColour  = util::paletteDarkNeutralColour();
            }
        } // colours

        // --------------------------------------------------------------------

        // init visualizer
        vis::MyVisPtr vptr( new pcl::visualization::PCLVisualizer(title) );
        vptr->setBackgroundColor( bg_colour(0), bg_colour(1), bg_colour(2) );

        // --------------------------------------------------------------------

        // create PCL pointcloud from input "points"
        MyPCLCloud::Ptr                      cloud  ( new MyPCLCloud );
        pcl::PointCloud<pcl::Normal>::Ptr normals( new pcl::PointCloud<pcl::Normal> );
        {
            cloud->reserve( points.size() );

            // Allocate PCL point
            MyPCLPoint pnt;

            // for each point in input
            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                const int                     point_gid      = points[pid].getTag(PointPrimitiveT::TAGS::GID); // which patch is the point assigned to
                GidLidLid1Map::const_iterator gidLidIterator = gid2lidLid1.find( point_gid );            // which primitive is in this patch ( only one...:( )
                LidLid1                       lidLid1(-1,-1);                                            // exact indices of primitive
                Colour                        point_colour   = unusedPointColour;                        // point colour

                // get primitive of point
                if ( gidLidIterator != gid2lidLid1.end() )
                    lidLid1 = gidLidIterator->second;

                // retrieve point colour
                if ( gidLidIterator != gid2lidLid1.end() )
                {
                    const int mid = primitives[ lidLid1.first ][ lidLid1.second ].getTag( primColourTag );
                    if ( id2ColId.find( mid ) != id2ColId.end() )
                    {
                        point_colour = pointColours[id2ColId[mid]];
                    }
                }

                // skip, if gid filtering is on, and this gid is not in the filter
                if ( filter_gids && (filter_gids->find( point_gid ) == filter_gids->end()) )
                    continue;

                // Set PCL point position
                if ( draw_mode & DRAW_MODE::REPROJECT ) // reprojected draw
                {
                    // cache point position
                    Position const& pos = points[pid].template pos();

                    // check, if point is assigned to existing patch
                    if ( gidLidIterator != gid2lidLid1.end() )
                        pnt.getVector3fMap() = primitives[ lidLid1.first ][ lidLid1.second ].projectPoint( pos );
                    // show point unprojected, NOT HIDE_UNASSIGNED_PTS
                    else if ( !(draw_mode & HIDE_UNASSIGNED_PTS) )
                        pnt.getVector3fMap() = pos;
                }
                else                            // regular draw
                {
                    pnt.getVector3fMap() = points[pid].template pos();
                }

                // Set PCL point colour
                pnt.r = point_colour( 0 );
                pnt.g = point_colour( 1 );
                pnt.b = point_colour( 2 );

                // Add PCL point
                cloud->push_back( pnt );

                // Add point normal, will be used, if needed
                if ( show_normals )
                {
                    pcl::Normal normal;
                    normal.normal_x = points[pid].template dir()(0);
                    normal.normal_y = points[pid].template dir()(1);
                    normal.normal_z = points[pid].template dir()(2);
                    normals->push_back( normal );
                } // ...show normals
            } //...for points
        } //...add PCL pointcloud

        // --------------------------------------------------------------------

        if ( !hide_points )
            vptr->addPointCloud( cloud, "cloud", 0 );

        // --------------------------------------------------------------------

        if ( show_normals )
        {
            vptr->addPointCloudNormals<MyPCLPoint,pcl::Normal>( /*  point_cloud: */ cloud
                                                              , /* normal_cloud: */ normals
                                                              , /*        level: */ show_normals // show every level-th normal
                                                              , /*        scale: */ 0.02f
                                                              , /*   cloud_name: */ "normal_cloud"
                                                              , /*  viewport_id: */ 0 );
        }

        // --------------------------------------------------------------------

        // show point ids (takes quite long time)
        if ( show_pids )
        {
            for ( int pid = 0; pid != points.size(); ++pid )
            {
                char pname[255],ptext[255];
                sprintf( pname, "p%d", pid );
                sprintf( ptext, "%d", points[pid].getTag(PointPrimitiveT::TAGS::GID) );
                vptr->addText3D( ptext, cloud->at(pid), 0.005, cloud->at(pid).r/255.f, cloud->at(pid).g/255.f, cloud->at(pid).b/255.f, pname, 0 );
            }
        }

        // --------------------------------------------------------------------

        // set point size
        if ( !hide_points )
            vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, point_size, "cloud", 0 );

        // --------------------------------------------------------------------

        // count populations
        GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        processing::getPopulations( populations, points );

        // problem
#if 0
        typedef double OptScalar; // Mosek, and Bonmin uses double internally, so that's what we have to do...
        typedef qcqpcpp::OptProblem<OptScalar> OptProblemT;
        OptProblemT *p_problem = new qcqpcpp::OptProblem<OptScalar>();
        if ( !problemPath.empty() )
        {
            int err = EXIT_SUCCESS;
            err += p_problem->read( problemPath );
            if ( EXIT_SUCCESS != err )
            {
                std::cerr << "[" << __func__ << "]: " << "Could not read problem, exiting" << std::endl;
                exit(0);
            }
        }
#endif

        // spatial significance cache variable
        Eigen::Matrix<_Scalar,Eigen::Dynamic,1> area( 1, 1 );

        // --------------------------------------------------------------------

        pcl::PolygonMesh hull_mesh_accum, plane_mesh;
        MyPCLCloud plane_mesh_cloud;                   // cloud to save points, and then add back to mesh in the end

        // --------------------------------------------------------------------

        // draw primitives
        if ( !(draw_mode & DRAW_MODE::HIDE_PRIMITIVES) )
        {
            for ( size_t lid = 0; lid != primitives.size(); ++lid )
                for ( size_t lid1 = 0; lid1 != primitives[lid].size(); ++lid1 )
                {
                    // caching
                    PrimitiveT const& prim = primitives[lid][lid1];
                    const int gid     = prim.getTag( PrimitiveT::TAGS::GID     );
                    const int dir_gid = prim.getTag( PrimitiveT::TAGS::DIR_GID );

                    // status filtering
                    if ( filter_status && (*filter_status).find(primitives[lid][lid1].getTag(PrimitiveT::TAGS::STATUS)) == (*filter_status).end() )
                        continue;

                    // WTF? // TODO remove
                    //if (primitives[lid][lid1].getTag(PrimitiveT::TAGS::STATUS) == PrimitiveT::STATUS_VALUES::SMALL)
                    //    continue;

                    // GID filtering
                    if ( filter_gids && (filter_gids->find(gid) == filter_gids->end()) )
                    {
                        continue;
                    }

                    // DIR_GID filtering
                    if ( filter_dids && (filter_dids->find(dir_gid) == filter_dids->end()) )
                    {
                        continue;
                    }

                    // unique primitive name
                    char line_name[64];
                    sprintf( line_name, "line_%04lu_%04lu", lid, lid1 );

                    // cache colour
                    const int colour_id = prim.getTag( primColourTag );
                    Colour prim_colour = primColours[ id2ColId[colour_id] ] / 255.;
                    //std::cout << "reading colour " << prim_colour.transpose() << " as id2ColId[" << primColourTag << "] = "
                    //          << id2ColId[primColourTag] << std::endl;

                    // use assignments: if use tags, collect GID tagged point indices
                    std::vector<PidT> indices;
                    if ( use_tags )
                    {
                        for ( PidT pid = 0; pid != points.size(); ++pid )
                            if ( points[pid].getTag(PointPrimitiveT::TAGS::GID) == gid )
                                indices.push_back( pid );

                        // don't show unpopulated primitives
                        if ( skip_empty && !indices.size() ) continue;

                        if ( use_tags == 1 ) // mode2 means no ellipses
                            drawEllipse( vptr, cloud, indices, scale, gid, prim_colour );

                    } //...if use_tags

                    // DRAW
                    PrimitiveT::template draw<PointPrimitiveT>( /*   primitive: */ primitives[lid][lid1]
                                                              , /*      points: */ points
                                                              , /*   threshold: */ scale
                                                              , /*     indices: */ use_tags ? &indices : NULL
                                                              , /*      viewer: */ vptr
                                                              , /*   unique_id: */ line_name
                                                              , /*      colour: */ prim_colour(0), prim_colour(1), prim_colour(2)
                                                              , /* viewport_id: */ 0
                                                              , /*     stretch: */ stretch /* = 1.2 */
                                                              , /*       qhull: */ old_draw_mode /* = 1, classic, axis aligned */
                                                              , /*       alpha: */ hull_alpha
                                                              );
                    // commented on 29/10/2014 by Aron, reason: should be in lineprimitive::draw
                    // vptr->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 4.0, line_name, 0 );

                    // store line size in "area"
                    if ( show_pop || show_ids )
                    {
                        primitives[lid][lid1].getSpatialSignificance( /* [out] sqrt(max(eigval)): */ area
                                                                    , /* [in]             points: */ points
                                                                    , /* [in]              scale: */ scale
                                                                    , /* [in]            indices: */ &(populations[gid]) );
                    }

                    // show patch size (spatial significance)
                    if ( show_pop && !lid1 ) // only once per cluster
                    {
                        char popstr[255];
                        sprintf( popstr, "%2.4f", area(0) );
                        vptr->addText3D( popstr, pclutil::asPointXYZ( primitives[lid][lid1].template pos() )
                                         , text_size, prim_colour(0), prim_colour(1), prim_colour(2), line_name + std::string("_pop"), 0 );
                    }

                    // show GID and DIR_GID
                    if ( show_ids )
                    {
                        char gid_name[255],gid_text[255];
                        sprintf( gid_name, "primgid%d_%d", gid, dir_gid  );
                        //sprintf( gid_text, "(%2.4f),%d,%d,%d", area(0), gid, dir_gid, primitives[lid][lid1].getTag(PrimitiveT::TAGS::STATUS) );
                        sprintf( gid_text, "(%2.4f),%d,%d,%.1f", area(0), gid, dir_gid, primitives[lid][lid1].getTag(PrimitiveT::TAGS::GEN_ANGLE) * deg_multiplier );
                        Position const& pos = primitives[lid][lid1].template pos();
                        vptr->addText3D( gid_text
                                       , pclutil::asPointXYZ( pos )
                                       , text_size //0.015
                                       , prim_colour(0)/2., prim_colour(1)/2., prim_colour(2)/2., gid_name, 0 );
                        if ( show_pop )
                            vptr->removeText3D( line_name + std::string("_pop"), 0 );
                    }

                    // draw connections
                    if ( angles || show_spatial )
                    {
                        // check for angles
                        if ( show_spatial && !angles )
                        {
                            std::cerr << "need angles to show distances" << std::endl;
                            throw new std::runtime_error("need angles");
                        }

                        typename PrimitiveT::ExtentsT extents0, extents1;
                        SpatialSqrtPrimitivePrimitiveEnergyFunctor<MyFinitePrimitiveToFinitePrimitiveCompatFunctor<PrimitiveT>,PointContainerT,Scalar,PrimitiveT>
                                distFunctor( *angles, points, scale );
                        distFunctor.setDirIdBias( 0 );
                        distFunctor.setSpatialWeightCoeff( 20 );
                        distFunctor.setTruncAngle( 0.3 );
                        distFunctor.setUseAngleGen( 1 );

                        if ( populations[gid].size() > pop_limit )
                        {
                            for ( size_t lid2 = lid; lid2 != primitives.size(); ++lid2 )
                            {
                                const int gid1 = primitives[lid2].at(0).getTag(PrimitiveT::TAGS::GID);
                                if ( populations[ gid1 ].size() < pop_limit )
                                    continue;

                                for ( size_t lid3 = lid1; lid3 < primitives[lid2].size(); ++lid3 )
                                {
                                    if ( (lid == lid2) && (lid1 == lid3) ) continue;
                                    PrimitiveT const& prim1 = primitives[lid2][lid3];
#if 1
                                    if ( show_spatial )
                                    {
                                        prim.template getExtent<PointPrimitiveT>( extents0
                                                                                , points
                                                                                , scale
                                                                                , &(populations[gid]) );
                                        prim1.template getExtent<PointPrimitiveT>( extents1
                                                                                , points
                                                                                , scale
                                                                                , &(populations[gid1]) );

                                        _Scalar idealAngle, spatW;
                                        _Scalar dist = distFunctor.eval( prim, extents0, prim1, extents1, *angles, &idealAngle, &spatW );
                                        if ( !(dist > _Scalar(0.)) ) continue;

                                        {
                                            char name[255];
                                            sprintf( name, "dist_l%lu%lu_l%lu%lu", lid, lid1, lid2, lid3 );
                                            vptr->addLine( pclutil::asPointXYZ( prim.pos() )
                                                         , pclutil::asPointXYZ( prim1.pos() )
                                                         , gray(0), gray(1), gray(2), name, 0 );
                                            vptr->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, name, 0 );

                                            pcl::PointXYZ line_center;
                                            line_center.getVector3fMap() = (prim.pos() + prim1.pos()) / _Scalar(2.);

                                            {
                                                char dist_str[255];
                                                if ( spatW > _Scalar(0.) )
                                                    sprintf( dist_str, "%.4f(%.0f,%.3f)", dist, idealAngle * deg_multiplier, spatW );
                                                else
                                                    sprintf( dist_str, "%.4f(%.0f)", dist, idealAngle * deg_multiplier );
                                                vptr->addText3D( dist_str
                                                                 , line_center, text_size, gray(0)+.1, gray(1)+.1, gray(2)+.1
                                                                 , name + std::string("_dist"), 0 );
                                            }
                                        }

                                        continue;
                                    }
#endif

                                    _Scalar angle = MyPrimitivePrimitiveAngleFunctor::eval( primitives[lid][lid1], primitives[lid2][lid3], *angles );
                                    if ( (angle < angle_limit) )
                                    {
                                        char name[255];
                                        sprintf( name, "conn_l%lu%lu_l%lu%lu", lid, lid1, lid2, lid3 );
                                        vptr->addLine( pclutil::asPointXYZ( primitives[lid][lid1].pos() )
                                                     , pclutil::asPointXYZ( primitives[lid2][lid3].pos() )
                                                     , gray(0), gray(1), gray(2), name, 0 );
                                        vptr->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, name, 0 );

                                        pcl::PointXYZ line_center;
                                        line_center.getVector3fMap() = (primitives[lid][lid1].pos() + primitives[lid2][lid3].pos()) / _Scalar(2.);

                                        if ( print_perf_angles )
                                        {
                                            char ang_str[255];
                                            sprintf( ang_str, "%.2f°", angle * deg_multiplier );
                                            vptr->addText3D( ang_str
                                                             , line_center, text_size, gray(0)+.1, gray(1)+.1, gray(2)+.1
                                                             , name + std::string("_ang"), 0 );
                                        }
                                    } //... if angle close enough
                                } //...for lid3
                            } //...for lid2
                        } //...if populations
                    } //...if angles

                    // red lines for same group id
                    if ( angles && !show_spatial )
                    {
                        for ( size_t lid2 = lid; lid2 != primitives.size(); ++lid2 )
                        {
                            for ( size_t lid3 = lid1; lid3 < primitives[lid2].size(); ++lid3 )
                            {
                                if ( (lid == lid2) && (lid1 == lid3) ) continue;

                                if ( primitives[lid2][lid3].getTag(PrimitiveT::TAGS::DIR_GID) == primitives[lid][lid1].getTag(PrimitiveT::TAGS::DIR_GID) )
                                {
                                    char name[255];
                                    sprintf( name, "same_l%lu%lu_l%lu%lu", lid, lid1, lid2, lid3 );
                                    vptr->addLine( pclutil::asPointXYZ( primitives[lid][lid1].pos() )
                                                   , pclutil::asPointXYZ( primitives[lid2][lid3].pos() )
                                                   , 1., .0, .0, name, 0 );
                                    vptr->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, name, 0 );
                                }
                            }
                        }
                    } //...red lines for same ids

                    // Polygon export
                    if ( save_poly && (PrimitiveT::EmbedSpaceDim == 3) )
                    {
                        // (1) QHULL
                        if ( draw_mode & DRAW_MODE::QHULL )
                        {
                            MyPCLCloud hull_cloud;
                            pcl::PolygonMesh mesh;
                            if ( PrimitiveT::getHull( hull_cloud, primitives[lid][lid1], points, &populations[gid], hull_alpha, &mesh ) )
                            {
                                std::cout << "mesh.cloud has " << mesh.cloud.width << " x " << mesh.cloud.height << " points";
                                const int pid_offs = hull_mesh_accum.cloud.width;
                                pcl::concatenatePointCloud( hull_mesh_accum.cloud, mesh.cloud, hull_mesh_accum.cloud );
                                hull_mesh_accum.polygons.resize( hull_mesh_accum.polygons.size() + 1 );
                                for ( int i = 0; i != mesh.polygons.size(); ++i )
                                {

                                    std::cout<<"mesh.polygons["<<i<<"].vertices:";
                                    for(size_t vi=0;vi!=mesh.polygons[i].vertices.size();++vi)
                                    {
                                        hull_mesh_accum.polygons.back().vertices.push_back( pid_offs + mesh.polygons[i].vertices[vi] );
                                        std::cout << mesh.polygons[i].vertices[vi] << "(" << pid_offs + mesh.polygons[i].vertices[vi] <<"), "; fflush(stdout);
                                    }
                                    std::cout << "\n";
                                }
                            }
                        } //...(1) QHull

                        // (2) Planes
                        {
                            std::vector<Position> minMax;
                            int err2 = primitives[lid][lid1].template getExtent<PointPrimitiveT>( /*             extent: */ minMax
                                                                                                , /*             points: */ points
                                                                                                , /*              scale: */ scale
                                                                                                , /*            indices: */ &(populations[gid])
                                                                                                , /* force_axis_aligned: */ draw_mode & DRAW_MODE::AXIS_ALIGNED );
                            // if extent exists
                            if ( err2 == EXIT_SUCCESS )
                            {
                                plane_mesh.polygons.resize( plane_mesh.polygons.size() + 1 );
                                const int pid_offs = plane_mesh_cloud.size();
                                for ( int i = 0; i != minMax.size(); ++i )
                                {
                                    MyPCLPoint pnt;
                                    pnt.getVector3fMap() = minMax[i];
                                    plane_mesh_cloud.push_back( pnt );
                                    plane_mesh      .polygons.back().vertices.push_back( pid_offs + i );
                                }
                            }
                        } //...(2) Planes
                    } //...polygon export
                } //...lid1
        } //...if (!draw_mode & HIDE_PRIMITIVES)

        // --------------------------------------------------------------------

        if ( !no_scale_sphere )
        {
            MyPCLPoint min_pt, max_pt;
            pcl::getMinMax3D( *cloud, min_pt, max_pt );

            vptr->addSphere( pcl::PointXYZ(0,0,0), scale, "scale_sphere", 0 );
        }

        // --------------------------------------------------------------------

        if ( save_poly && (PrimitiveT::EmbedSpaceDim == 3) )
        {
            if ( draw_mode & DRAW_MODE::QHULL )
            {
                std::cout << "mesh_accum.size: " << hull_mesh_accum.cloud.width << ", and " << hull_mesh_accum.polygons.size() << " polygons" << std::endl;
                pcl::io::savePLYFile( "mesh.ply", hull_mesh_accum );
            }
            pcl::toPCLPointCloud2( plane_mesh_cloud, plane_mesh.cloud );
            std::cout << "plane_mesh.size: " << plane_mesh.cloud.width << ", and " << plane_mesh.polygons.size() << " polygons" << std::endl;
            pcl::io::savePLYFile( "plane_mesh.ply", plane_mesh );
        }

        if ( save_hough )
        {
            std::string     houghFilePath( "hough.csv"   );
            std::ofstream   houghFile    ( houghFilePath );
            if ( !houghFile.is_open() )
            {
                std::cerr << "[" << __func__ << "]: " << "could not open " << houghFilePath << std::endl;
                throw new std::runtime_error("");
            }

            for ( size_t lid = 0; lid != primitives.size(); ++lid )
                for ( size_t lid1 = 0; lid1 < primitives[lid].size(); ++lid1 )
                {
                    PrimitiveT const& prim = primitives[lid][lid1];

                    Eigen::Matrix<Scalar,3,1> dir = prim.dir();
                    houghFile << angleInRad( dir, Eigen::Matrix<Scalar,3,1>::UnitX().eval() ) << ","
                              << prim.getDistance( Eigen::Matrix<Scalar,3,1>::Zero() )
                              << std::endl;
                }

            houghFile.close();
        }

        // --------------------------------------------------------------------

        if ( spin )
            vptr->spin();
        else
            vptr->spinOnce();

        return vptr;
    } // ...Visualizer::show()

    template <class PrimitiveContainerT, class PointContainerT>
    template <typename _Scalar> int
    Visualizer<PrimitiveContainerT,PointContainerT>::drawEllipse( vis::MyVisPtr                             vptr
                                                                , MyPCLCloud::Ptr                              cloud
                                                                , std::vector<PidT>                  const& indices
                                                                , _Scalar                            const  scale
                                                                , int                                const  prim_tag
                                                                , Eigen::Matrix<_Scalar,3,1>         const& prim_colour
                                                                )
    {
        typedef typename PrimitiveContainerT::value_type::value_type PrimitiveT;
        typedef typename PointContainerT::value_type PointT;

        pcl::PCA<MyPCLPoint> pca;
        pca.setInputCloud( cloud );
        pcl::PointIndices::Ptr indices_ptr( new pcl::PointIndices() );
        indices_ptr->indices.resize( indices.size() );
        std::copy( indices.begin(), indices.end(), indices_ptr->indices.begin() );
        //indices_ptr->indices = indices;
        pca.setIndices( indices_ptr );

        const _Scalar min_dim1 = scale * 0.25f;
        MyPCLCloud::Ptr poly( new MyPCLCloud );
        Eigen::Matrix<_Scalar,3,1> centroid( Eigen::Matrix<_Scalar,3,1>::Zero() );
        Eigen::Matrix<_Scalar,2,1> dims; dims << 0, 0;
        Eigen::Matrix<_Scalar,3,3> eigen_vectors( Eigen::Matrix<_Scalar,3,3>::Zero() );
        if ( indices.size() > 2 )
        {
            centroid = pca.getMean().head<3>();
            eigen_vectors = pca.getEigenVectors();

            for ( int pid_id = 0; pid_id != indices.size(); ++pid_id )
            {
                const int pid = indices[ pid_id ];
                MyPCLPoint pnt;
                pca.project( cloud->at(pid), pnt );
                if ( std::abs(pnt.x) > dims(0) ) dims(0) = std::abs(pnt.x);
                if ( std::abs(pnt.y) > dims(1) ) dims(1) = std::abs(pnt.y);
            }
            if ( dims(1) < min_dim1 )
                dims(1) = min_dim1;
        }
        else if ( indices.size() == 2 ) // line
        {
            centroid = (cloud->at(indices[1]).getVector3fMap() + cloud->at(indices[0]).getVector3fMap()) / _Scalar(2);
            eigen_vectors.col(0) = (cloud->at(indices[1]).getVector3fMap() - centroid).normalized();
            eigen_vectors.col(1) = eigen_vectors.col(0).cross( Eigen::Matrix<_Scalar,3,1>::UnitZ() );
            dims(0) = std::max( (cloud->at(indices[1]).getVector3fMap() - centroid).norm(), min_dim1 );
            dims(1) = std::min( dims(0)/2.f, min_dim1 );
        }

        char poly_name[255];
        static int sphere_id = 0;
        sprintf( poly_name, "poly%d_%d", prim_tag, sphere_id++ );

        if ( indices.size() > 1 )
        {
            MyPCLPoint tmp;
            for ( int dir = -1; dir <= 1; dir += 2 )
                for ( int dim = 0; dim != 2; ++dim )
                {
                    tmp.getVector3fMap() = centroid + (_Scalar(dir) * eigen_vectors.col(dim) * dims(dim)) * 1.1f;
                    poly->push_back( tmp );
                }

                vptr->addPolygon<MyPCLPoint>( poly
                                            , prim_colour(0), prim_colour(1), prim_colour(2)
                                            , poly_name
                                            , 0 );
//            vptr->addText3D( poly_name, poly->at(0), 0.01,  prim_colour(0), prim_colour(1), prim_colour(2)
//            , std::string(poly_name)+"text", 0 );
        }
        else if ( indices.size() == 1 )
        {
            pcl::ModelCoefficients circle_coeffs;
            circle_coeffs.values.resize(3);
            circle_coeffs.values[0] = cloud->at(indices[0]).x;
            circle_coeffs.values[1] = cloud->at(indices[0]).y;
            circle_coeffs.values[2] = min_dim1;
            vptr->addCircle( circle_coeffs, poly_name, 0 );
        }

        return 0;
    } // ... Visulazier::drawEllipse
} // ns gf2

#endif // __GF2_VISUALIZER_HPP__
