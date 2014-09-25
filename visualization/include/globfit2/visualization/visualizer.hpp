#ifndef __GF2_VISUALIZER_HPP__
#define __GF2_VISUALIZER_HPP__

#if GF2_USE_PCL
#   include "pcl/point_types.h"
#   include "pcl/visualization/pcl_visualizer.h"
#endif

#include "globfit2/visualization/visualization.h" // MyVisPtr
#include "globfit2/processing/util.hpp"           // getPopulations()

namespace GF2 {
    //! \brief Visualizer class to show points, primitives and their relations
    //! \tparam PrimitiveContainerT     Patch-wise grouped primitive container type. Concept: vector< vector< GF2::LinePrimitive2> >.
    //! \tparam PointContainerT         Cloud of possibly oriented 3D points. Concept: vector< GF2::PointPrimitive >.
    template <class PrimitiveContainerT, class PointContainerT>
    class Visualizer
    {
            typedef pcl::PointXYZRGB         MyPoint;
            typedef pcl::PointCloud<MyPoint> MyCloud;
        public:
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
             *  \return             The visualizer for further display and manipulation
             */
            template <typename _Scalar> static inline vis::MyVisPtr
            show( PrimitiveContainerT  const& primitives
                , PointContainerT      const& points
                , _Scalar              const  scale
                , Eigen::Vector3f      const& colour                = (Eigen::Vector3f() << 0.f, 0.f, 1.f).finished()
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
                , std::set<int>        const* filter_status         = NULL
                , _Scalar              const  stretch               = 1.
                , int                  const  draw_mode             = 0
                );

            //! \brief Shows a polygon that approximates the bounding ellipse of a cluster
            template <typename _Scalar> static inline int
            drawEllipse( vis::MyVisPtr                             vptr
                       , MyCloud::Ptr                              cloud
                       , std::vector<int>                   const& indices
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
                                                           , Eigen::Vector3f      const& colour              /* = {0,0,1} */
                                                           , bool                 const  spin                /* = true */
                                                           , std::vector<_Scalar> const* angles              /* = NULL */
                                                           , bool                 const  show_ids            /* = false */
                                                           , char                 const  use_tags            /* = false */
                                                           , int                  const  pop_limit           /* = 10 */
                                                           , std::string          const& title               /* = "" */
                                                           , bool                 const  show_pids           /* = false */
                                                           , int                  const  show_normals        /* = 0 */
                                                           , bool                 const  show_pop            /* = false */
                                                           , _Scalar              const  perfect_angle_limit /* = 0872664625 = 5deg*/
                                                           , bool                 const  print_perf_angles   /* = false */
                                                           , bool                 const  dir_colours         /* = false */
                                                           , bool                 const  hide_points         /* = false */
                                                           , std::set<int>        const* filter_gids         /* = NULL */
                                                           , std::set<int>        const* filter_status       /* = NULL */
                                                           , _Scalar              const  stretch             /* = 1. */
                                                           , int                  const  draw_mode           /* = 0 */
                                                           )
    {
#if 1
        typedef typename PrimitiveContainerT::value_type::value_type PrimitiveT;
        typedef typename PointContainerT::value_type PointPrimitiveT;

        const Eigen::Matrix<double,3,1> gray          ( (Eigen::Matrix<double,3,1>()<<.6,.6,.5).finished() );
        const _Scalar                   deg_multiplier( _Scalar(180.) / M_PI );

        std::cout << "[" << __func__ << "]: " << "scale: " << scale
                  << std::endl;

        const int primColourTag  = dir_colours ? PrimitiveT::DIR_GID      : PrimitiveT::GID;
        const int pointColourTag = dir_colours ? PointPrimitiveT::DIR_GID : PointPrimitiveT::GID;

        // calc groups
        int max_gid = 0, nPrimitives = 0, max_dir_gid = 0, nbColour = 0;
        std::map< int, std::pair<int,int> > gid2lidLid1;
        std::map< int, int >  id2ColId;
        for ( size_t pid = 0; pid != points.size(); ++pid )
            max_gid = std::max( max_gid, points[pid].getTag(PointPrimitiveT::GID) );
        for ( size_t lid = 0; lid != primitives.size(); ++lid )
            for ( size_t lid1 = 0; lid1 != primitives[lid].size(); ++lid1 )
            {
                if (primitives[lid][lid1].getTag(PrimitiveT::STATUS) != PrimitiveT::SMALL){

                    max_gid     = std::max( max_gid    , primitives[lid][lid1].getTag(PrimitiveT::GID    ) );
                    max_dir_gid = std::max( max_dir_gid, primitives[lid][lid1].getTag(PrimitiveT::DIR_GID) );

                    // this id has not been referenced to a color previously
                    if (id2ColId.find(primitives[lid][lid1].getTag(primColourTag)) == id2ColId.end()){
                        id2ColId  [ primitives[lid][lid1].getTag(primColourTag) ] = nbColour;
                        nbColour++;
                    }
                    gid2lidLid1[primitives[lid][lid1].getTag(PrimitiveT::GID)] = std::pair<int,int>(lid,lid1);

                    ++nPrimitives;
                }
            }
        std::cout << "[" << __func__ << "]: "
                  << "points: "         << points.size()
                  << ", primitives: "   << nPrimitives
                  << ", max_gid: "      << max_gid
                  << ", max_dir_gid: "  << max_dir_gid
                  << ", pop-limit: "    << pop_limit
                  << ", filter-gids.size(): " << (filter_gids ? filter_gids->size() : 0)
                  << std::endl;


        // Check if we can use color palette
        bool usePalette =  id2ColId.size() < util::paletteLightColoursEigen().size();
        std::vector<Eigen::Vector3f> pointColours, primColours;
        Eigen::Vector3f unusedPointColour, unusedPrimColour;

        if ( usePalette )
        {
            cout << "SWitching to nicer color palette" << endl;

            pointColours = util::paletteLightColoursEigen();
            primColours  = util::paletteDarkColoursEigen();

            unusedPointColour = util::paletteLightNeutralColour();
            unusedPrimColour  = util::paletteDarkNeutralColour();

        }
        else
        {
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
            unusedPointColour = unusedPrimColour = Eigen::Vector3f::Zero();
        }

        //pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
        vis::MyVisPtr vptr( new pcl::visualization::PCLVisualizer(title) );
        vptr->setBackgroundColor( 1, 1, 1 );
        MyCloud::Ptr cloud( new MyCloud );
        pcl::PointCloud<pcl::Normal>::Ptr normals( new pcl::PointCloud<pcl::Normal> );
        {
            cloud->reserve( points.size() );
            for ( size_t pid = 0; pid != points.size(); ++pid )
            {

                Eigen::Vector3f colour = unusedPointColour;

                if ( gid2lidLid1.find(points[pid].getTag(PointPrimitiveT::GID)) != gid2lidLid1.end() )
                {
                    std::pair<int,int> lid =  gid2lidLid1[points[pid].getTag(PointPrimitiveT::GID)];
                    std::cout << "pid" << pid << ".lid=" << lid.first << ", " << lid.second << ", primgid: " << primitives[lid.first][lid.second].getTag(PrimitiveT::GID)
                            << ", primdid: " << primitives[lid.first][lid.second].getTag(PrimitiveT::DIR_GID)
                            << ", pgid: " << points[pid].getTag(PointPrimitiveT::GID) << std::endl;
                    const int mid = primitives[lid.first][lid.second].getTag(primColourTag);
                    std::cout << ", mid= " << mid << ", for primcolourtag: " << primColourTag;
                    if ( id2ColId.find(mid) != id2ColId.end() )
                    {
                        if ( primitives[lid.first][lid.second].getTag(PrimitiveT::GID) == 15 )
                            std::cout << "here";
                        colour = pointColours[id2ColId[mid]];
                        std::cout << ", found colour " << id2ColId[mid] << ": "  << colour.transpose() << std::endl;
                    }
                    else
                        std::cout << ", didn't fing colour, so " << colour.transpose() << std::endl;
                }
                // skip, if filtering is on, and this gid is not in the filter
                //if ( filter_gids && (filter_gids->find(primitives[points[pid]].getTag(PrimitiveT::GID)) == filter_gids->end()) )
                //{
                //    continue;
                //}

                MyPoint pnt;
                pnt.x = ((Eigen::Matrix<_Scalar,3,1> )points[pid])(0); // convert PointPrimitive to Eigen::Matrix, and get (0)
                pnt.y = ((Eigen::Matrix<_Scalar,3,1> )points[pid])(1);
                pnt.z = ((Eigen::Matrix<_Scalar,3,1> )points[pid])(2);
                pnt.r = colour(0);
                pnt.g = colour(1);
                pnt.b = colour(2);
                cloud->push_back( pnt );

                pcl::Normal normal;
                normal.normal_x = points[pid].template dir()(0);
                normal.normal_y = points[pid].template dir()(1);
                normal.normal_z = points[pid].template dir()(2);
                normals->push_back( normal );
            }
        }
        if (! hide_points)
            vptr->addPointCloud( cloud, "cloud", 0 );
        if ( show_normals )
        {
            vptr->addPointCloudNormals<MyPoint,pcl::Normal>( cloud
                                                             , normals
                                                             , /* level: */ show_normals
                                                             , /* scale: */ 0.02f
                                                             , /* cloud_name: */ "normal_cloud"
                                                             , /* vid: */ 0 );
        }
        if ( show_pids )
        {
            for ( int pid = 0; pid != points.size(); ++pid )
            {
                char pname[255],ptext[255];
                sprintf( pname, "p%d", pid );
                sprintf( ptext, "%d", points[pid].getTag(PointPrimitiveT::GID) );
                vptr->addText3D( ptext, cloud->at(pid), 0.005, cloud->at(pid).r/255.f, cloud->at(pid).g/255.f, cloud->at(pid).b/255.f, pname, 0 );
            }
        }
        vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6.0, "cloud", 0 );

        // count populations
        GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        processing::getPopulations( populations, points );

        Eigen::Matrix<_Scalar,Eigen::Dynamic,1> area(1,1);

        for ( size_t lid = 0; lid != primitives.size(); ++lid )
            for ( size_t lid1 = 0; lid1 != primitives[lid].size(); ++lid1 )
            {
                if ( filter_status && (*filter_status).find(primitives[lid][lid1].getTag(PrimitiveT::STATUS)) == (*filter_status).end() )
                    continue;

                char line_name[64];
                sprintf( line_name, "line_%04lu_%04lu", lid, lid1 );
                const int gid     = primitives[lid][lid1].getTag( PrimitiveT::GID     );
                const int dir_gid = primitives[lid][lid1].getTag( PrimitiveT::DIR_GID );

                if (primitives[lid][lid1].getTag(PrimitiveT::STATUS) == PrimitiveT::SMALL) continue;

                //if ( filter_gids && (filter_gids->find(gid) == filter_gids->end()) )
                //{
                //    continue;
                //}

//                if ( lid != gid )
//                { std::cout << "!!\tstarting " << line_name << ": " << lid << ", " << lid1 << ", " << gid << ", " << dir_gid << std::endl; fflush(stdout); }
//                else
//                { std::cout << "starting " << line_name << ": " << lid << ", " << lid1 << ", " << gid << ", " << dir_gid << std::endl; fflush(stdout); }

//                if ( strcmp(line_name,"line_2_11") == 0 )
//                { std::cout << "!!!\t\tequal " << line_name << ": " << lid << ", " << lid1 << ", " << gid << ", " << dir_gid << std::endl; fflush(stdout); }

                Eigen::Matrix<_Scalar,3,1> prim_colour = primColours[id2ColId[dir_colours ? dir_gid : gid]] / 255.;
                if ( gid == 15 )
                    std::cout << "getting colour " << "id2ColId[" << (dir_colours ? dir_gid : gid)
                              << "] = " << id2ColId[dir_colours ? dir_gid : gid] << " = " << prim_colour.transpose() << " for 15 " << std::endl;

                // if use tags, collect GID tagged point indices
                std::vector<int> indices;
                if ( use_tags )
                {
                    for ( int pid = 0; pid != points.size(); ++pid )
                        if ( points[pid].getTag(PointPrimitiveT::GID) == gid )
                            indices.push_back( pid );
                    if ( use_tags == 1 ) // mode2 means no ellipses
                        drawEllipse( vptr, cloud, indices, scale, gid, prim_colour );

                } //...if use_tags
                PrimitiveT::template draw<PointPrimitiveT>( /*   primitive: */ primitives[lid][lid1]
                                                          , /*      points: */ points
                                                          , /*   threshold: */ scale
                                                          , /*     indices: */ use_tags ? &indices : NULL
                                                          , /*      viewer: */ vptr
                                                          , /*   unique_id: */ line_name
                                                          , /*      colour: */ prim_colour(0), prim_colour(1), prim_colour(2)
                                                          , /* viewport_id: */ 0
                                                          , /*     stretch: */ stretch /* = 1.2 */
                                                          , /*       qhull: */ draw_mode /* = 1, classic, axis aligned */
                                                          );
                vptr->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 4.0, line_name, 0 );

                // store line size in "area"
                if ( show_pop || show_ids )
                {
                    primitives[lid][lid1].getSpatialSignificance( /* [out] sqrt(max(eigval)): */ area
                                                                , /* [in]             points: */ points
                                                                , /* [in]              scale: */ scale
                                                                , /* [in]            indices: */ &(populations[gid]) );
                }

                if ( show_pop && !lid1 ) // only once per cluster
                {
                    char popstr[255];
                    //sprintf( popstr, "%d", populations[gid] );
                    sprintf( popstr, "%2.4f", area(0) );
                    vptr->addText3D( popstr, pclutil::asPointXYZ( primitives[lid][lid1].template pos() )
                                     , 0.015, prim_colour(0), prim_colour(1), prim_colour(2), line_name + std::string("_pop"), 0 );
                }

                // show line gid and dir gid
                if ( show_ids )
                {
                    char gid_name[255],gid_text[255];
                    sprintf( gid_name, "primgid%d_%d", gid, dir_gid  );
                    sprintf( gid_text, "(%2.4f),%d,%d,%d", area(0), gid, dir_gid, primitives[lid][lid1].getTag(PrimitiveT::STATUS) );
                    Eigen::Matrix<_Scalar,3,1> pos = primitives[lid][lid1].template pos();// + (Eigen::Matrix<_Scalar,3,1>() << scale, scale, 0.).finished();
                    vptr->addText3D( gid_text
                                   , pclutil::asPointXYZ( pos )
                                   , 0.015 //0.02
                                   , prim_colour(0)/2., prim_colour(1)/2., prim_colour(2)/2., gid_name, 0 );
                    if ( show_pop )
                    {
                        vptr->removeText3D( line_name + std::string("_pop"), 0 );
                    }
                }

                // draw connections
                if ( angles )
                {
                    if ( populations[gid].size() > pop_limit )
                    {
                        for ( size_t lid2 = lid; lid2 != primitives.size(); ++lid2 )
                        {
                            if ( populations[ primitives[lid2].at(0).getTag(PrimitiveT::GID) ].size() < pop_limit )
                                continue;

                            for ( size_t lid3 = lid1; lid3 < primitives[lid2].size(); ++lid3 )
                            {
                                if ( (lid == lid2) && (lid1 == lid3) ) continue;

                                _Scalar angle = MyPrimitivePrimitiveAngleFunctor::eval( primitives[lid][lid1], primitives[lid2][lid3], *angles );
                                if ( (angle < perfect_angle_limit) )
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
                                                         , line_center, 0.015, gray(0)+.1, gray(1)+.1, gray(2)+.1, name + std::string("_ang"), 0 );
                                    }
                                } //... if angle close enough
                            } //...for lid3
                        } //...for lid2
                    } //...if populations
                } //...if angles

                // red lines for same group
                if ( angles )
                {
                    for ( size_t lid2 = lid; lid2 != primitives.size(); ++lid2 )
                    {
                        for ( size_t lid3 = lid1; lid3 < primitives[lid2].size(); ++lid3 )
                        {
                            if ( (lid == lid2) && (lid1 == lid3) ) continue;

                            if ( primitives[lid2][lid3].getTag(PrimitiveT::DIR_GID) == primitives[lid][lid1].getTag(PrimitiveT::DIR_GID) )
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
                }
            } // ... lid1

        MyPoint min_pt, max_pt;
        pcl::getMinMax3D( *cloud, min_pt, max_pt );
        vptr->addSphere( pcl::PointXYZ(0,0,0), scale, "scale_sphere", 0 );
        if ( spin )
            vptr->spin();
        else
            vptr->spinOnce();

        return vptr;
#endif
        return vis::MyVisPtr();
    } // ...Visualizer::show()

    template <class PrimitiveContainerT, class PointContainerT>
    template <typename _Scalar> int
    Visualizer<PrimitiveContainerT,PointContainerT>::drawEllipse( vis::MyVisPtr                             vptr
                                                                , MyCloud::Ptr                              cloud
                                                                , std::vector<int>                   const& indices
                                                                , _Scalar                            const  scale
                                                                , int                                const  prim_tag
                                                                , Eigen::Matrix<_Scalar,3,1>         const& prim_colour
                                                                )
    {
#if 1
        typedef typename PrimitiveContainerT::value_type::value_type PrimitiveT;
        typedef typename PointContainerT::value_type PointT;

        pcl::PCA<MyPoint> pca;
        pca.setInputCloud( cloud );
        pcl::PointIndices::Ptr indices_ptr( new pcl::PointIndices() );
        indices_ptr->indices = indices;
        pca.setIndices( indices_ptr );

        const _Scalar min_dim1 = scale * 0.25f;
        MyCloud::Ptr poly( new MyCloud );
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
                MyPoint pnt;
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
            MyPoint tmp;
            for ( int dir = -1; dir <= 1; dir += 2 )
                for ( int dim = 0; dim != 2; ++dim )
                {
                    tmp.getVector3fMap() = centroid + (_Scalar(dir) * eigen_vectors.col(dim) * dims(dim)) * 1.1f;
                    poly->push_back( tmp );
                }

                vptr->addPolygon<MyPoint>( poly
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

#endif
        return 0;
    } // ... Visulazier::drawEllipse
} // ns gf2

#endif // __GF2_VISUALIZER_HPP__
