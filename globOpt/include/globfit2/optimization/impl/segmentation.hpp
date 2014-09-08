#ifndef GF2_SEGMENTATION_HPP
#define GF2_SEGMENTATION_HPP

#include "globfit2/optimization/segmentation.h"
#include "globfit2/my_types.h"                  // PCLPointAllocator
#include <vector>

#if GF2_USE_PCL
#   include "pcl/point_types.h"
#   include <pcl/point_cloud.h>
#endif

#include "globfit2/util/containers.hpp" // add( map, gid, primitive), add( vector, gid, primitive )
#include "globfit2/parameters.h"        // CandidateGeneratorParams

namespace GF2 {

//! \param[in,out] points
template < class    _PointPrimitiveT
         , class    _PrimitiveT
         , typename _Scalar
         , class    _PointContainerT
         >
int
Segmentation::orientPoints( _PointContainerT          &points
                          , _Scalar             const  scale
                          , int                 const  nn_K   )
{
    typedef pcl::PointCloud<pcl::PointXYZ>        CloudXYZ;

    // to pcl cloud
    CloudXYZ::Ptr cloud( new CloudXYZ() );
    _PointPrimitiveT::template toCloud<CloudXYZ::Ptr, _PointContainerT, PCLPointAllocator<_PointPrimitiveT::Dim> >
            ( cloud, points );

    // (1) local fit lines from pcl cloud
    std::vector<_PrimitiveT> fit_lines;
    std::vector<int        > point_ids;
    {
        fitLocal( /* [out]       lines: */ fit_lines
               , /*            points: */ cloud
               , /*           indices: */ NULL
               , /*              nn_K: */ nn_K
               , /*         nn_radius: */ scale
               , /*       soft_radius: */ true
               , /* [out]     mapping: */ &point_ids ); // contains point id for fit_line

        // copy line direction into point
        for ( int pid_id = 0; pid_id != point_ids.size(); ++pid_id )
        {
            const int pid = point_ids[pid_id];
            points[pid].coeffs().template segment<3>(3) = fit_lines[pid_id].dir();
        }
    } // ... (1) local fit

    return EXIT_SUCCESS;
} //...Segmentation::orientPoints()

template < class _PrimitiveContainerT
         , class _PointContainerPtrT> int
Segmentation:: fitLocal( _PrimitiveContainerT        & lines
                       , _PointContainerPtrT    const  cloud
                       , std::vector<int>       const* indices
                       , int                    const  K
                       , float                  const  radius
                       , bool                   const  soft_radius
                       , std::vector<int>            * point_ids
                       )
{
    using std::vector;
    typedef typename _PrimitiveContainerT::value_type         TLine;
    typedef typename TLine::Scalar              Scalar;
    typedef typename _PointContainerPtrT::element_type   PointsT;

    if ( indices ) { std::cerr << __PRETTY_FUNCTION__ << "]: indices must be NULL, not implemented yet..." << std::endl; return EXIT_FAILURE; }

    // get neighbourhoods
    std::vector< std::vector<int   > > neighs;
    std::vector< std::vector<Scalar> > sqr_dists;
    smartgeometry::getNeighbourhoodIndices( /*   [out] neighbours: */ neighs
                                            , /* [in]  pointCloud: */ cloud
                                            , /* [in]     indices: */ indices
                                            , /* [out]  sqr_dists: */ &sqr_dists
                                            , /* [in]        nn_K: */ K              // 5
                                            , /* [in]      radius: */ radius         // 0.02f
                                            , /* [in] soft_radius: */ soft_radius    // true
                                            );

    // only use, if more then 2 data-points
    if ( std::count_if( neighs.begin(), neighs.end(), [] (vector<int> const& n1) { return n1.size() > 2; } ) < 2 )
    {
        std::cerr << "[" << __func__ << "]: " << "not enough to work with (<2)...change scale " << radius << std::endl;
        return EXIT_SUCCESS;
    }

    // every point proposes primitive[s] using its neighbourhood
    int skipped = 0;
    for ( size_t pid = 0; pid != neighs.size(); ++pid )
    {
        // can't fit a line to 0 or 1 points
        if ( neighs[pid].size() < 2 )
        {
            ++skipped;
            continue;
        }

        Eigen::Matrix<Scalar,TLine::Dim,1> line;
        int err = smartgeometry::geometry::fitLinearPrimitive<PointsT,Scalar,TLine::Dim>( /*           output: */ line
                                                                                          , /*         points: */ *cloud
                                                                                          , /*          scale: */ radius
                                                                                          , /*        indices: */ &(neighs[pid])
                                                                                          , /*    refit times: */ 2
                                                                                          , /* use input line: */ false
                                                                                          );
        if ( err == EXIT_SUCCESS )      lines.emplace_back( TLine(line) );
        if ( point_ids )
        {
            point_ids->emplace_back( pid );
        }
    }
    std::cout << "[" << __func__ << "]: "
              << skipped << "/" << neighs.size() << ": " << skipped / static_cast<float>(neighs.size()) * 100.f << "% of points did not produce primitives, so the primitive count is:"
              << lines.size() << " = " << lines.size() / static_cast<float>(neighs.size()) *100.f << "%" << std::endl;

    return EXIT_SUCCESS;
} // ...Segment::propose()

/*
 * \brief Groups unoriented points into oriented patches represented by a single primitive
 *                   (1) group to patches
 *                   (2) refit lines to patches
 * \tparam _PrimitiveContainerT Concept: std::map< int, std::vector<_PrimitiveT> >. Groups primitives by their GID.
 * \param[out] patches                   Ready to process patches, each of them with one direction only. ( { 0: [Primitive00] }, { 1: [Primitive11] }, ... )
 * \param[in/out] points                 Input points that get assigned to the patches by their GID tag set.
 * \param[in]  scale                     Spatial scale to use for fits.
 * \param[in]  angles                    Desired angles to use for groupings.
 * \param[in]  patchPatchDistanceFunctor #regionGrow() uses the thresholds encoded to group points. The evalSpatial() function is used to assign orphan points.
 * \param[in]  nn_K                      Number of nearest neighbour points looked for in #regionGrow().
 */
template < class       _PrimitiveT
         , typename    _Scalar
         , class       _PrimitiveContainerT
         , class       _PointContainerT
         , class       _PatchPatchDistanceFunctorT
         > int
Segmentation::patchify( _PrimitiveContainerT                   & patches
                      , _PointContainerT                       & points
                      , _Scalar                           const  scale
                      , std::vector<_Scalar>              const& angles
                      , _PatchPatchDistanceFunctorT       const& patchPatchDistanceFunctor
                      , int                               const  nn_K
                      )
{
    typedef segmentation::Patch<_Scalar,_PrimitiveT> PatchT;
    typedef std::vector< PatchT >                    PatchesT;

    typedef typename _PointContainerT::value_type PointPrimitiveT;

    // log
    std::cout << "[" << __func__ << "]: " << "PatchPatchDistance by " << patchPatchDistanceFunctor.toString() << std::endl;

    // (1) group
    PatchesT groups;
    {
        regionGrow<_PrimitiveT>
                  ( /* [in,out]  points/pointsWGIDTag: */ points
                  , /* [out]          groups_pointids: */ groups
                  , /* [in]                     scale: */ scale
                  , /* [in] patchPatchDistanceFunctor: */ patchPatchDistanceFunctor
                  , /* [in]              gid_tag_name: */ PointPrimitiveT::GID
                  , /* [in]                      nn_K: */ nn_K );
    } // ... (1) group

    // (2) Create PrimitiveContainer
    // Copy the representative direction of each patch in groups to an output patch with GID as it's linear index in groups.
    for ( size_t gid = 0; gid != groups.size(); ++gid )
    {
        containers::add( patches, gid, groups[gid].getRepresentative() )
                .setTag( _PrimitiveT::GID    , gid )
                .setTag( _PrimitiveT::DIR_GID, gid );
    }

    return EXIT_SUCCESS;
} // ...Segmentation::patchify()

/*  \brief                               Greedy region growing
 *  \tparam _PrimitiveContainerT         Concept: std::vector<\ref GF2::LinePrimitive2>
 *  \tparam _PointContainerT             Concept: std::vector<\ref GF2::PointPrimitive>
 *  \tparam _PointPatchDistanceFunctorT  Concept: \ref RepresentativeSqrPatchPatchDistanceFunctorT.
 *  \tparam _PatchesT                    Concept: vector< \ref segmentation::Patch <_Scalar,_PrimitiveT> >
 *  \tparam _PrimitiveT                  Concept: \ref GF2::LinePrimitive2
 *  \tparam _Scalar                      Concept: float
 *  \tparam _PointT                      Concept: \ref GF2::PointPrimitive
 *  \param[in,out] points                Input points to create patches from. The \p gid_tag_name field of the points will be set according to their patch assignment.
 *  \param[out] groups_arg               Holds the point-id-groups, that can then be refit to to get a patch location and direction. Concept: vector<\ref segmentation::Patch>.
 *  \param[in] scale                     Spatial extent of the input. In practice unused, since the \p patchPatchDistanceFunctor was constructed with it.
 *  \param[in] patchPatchDistanceFunctor Takes two patches, and decides, whether they are similar enough to be merged.
 *                                       In practice, takes a patch that is currently grown, and a temporary patch
 *                                       that only contains a neighbouring point, and decides. See in \ref RepresentativeSqrPatchPatchDistanceFunctorT.
 *  \param[in] gid_tag_name              The key value of GID in _PointT. Suggested to be: _PointT::GID.
 *  \param[in] nn_K                      Number of nearest neighbour points looked for.
 */
template < class       _PrimitiveT
         , class       _PointContainerT
         , class       _PatchPatchDistanceFunctorT
         , class       _PatchesT
         , typename    _Scalar
         , class       _PointPrimitiveT> int
Segmentation::regionGrow( _PointContainerT                       & points
                        , _PatchesT                              & groups_arg
                        , _Scalar                           const  /*scale*/
                        , _PatchPatchDistanceFunctorT       const& patchPatchDistanceFunctor
                        , int                               const  gid_tag_name
                        , int                               const  nn_K )
{
    std::cout << "[" << __func__ << "]: " << "running with " << patchPatchDistanceFunctor.toString() << std::endl;
    std::cout << "[" << __func__ << "]: " << "running at " << patchPatchDistanceFunctor.getSpatialThreshold() << " spatial threshold" << std::endl;
    std::cout << "[" << __func__ << "]: " << "running at " << patchPatchDistanceFunctor.getAngularThreshold() << " radius threshold" << std::endl;

    typedef typename _PatchesT::value_type   PatchT;
    typedef          std::vector<PatchT>     Patches;

    // create patches with a single point in them
    std::deque<int> starting_cands( points.size() );
    for ( int pid = 0; pid != points.size(); ++pid )
    {
        //const int pid = point_ids_arg ? (*point_ids_arg)[ pid_id ] : pid_id;
        if ( points[pid].template dir().norm() != _Scalar(1) )
            std::cerr << "unoriented point at pid " << pid << "? " << points[pid].template dir().transpose() << ", norm: " << points[pid].template dir().norm() << std::endl;
        starting_cands.push_back( pid );
    }

    // prebulid ann cloud
    pcl::PointCloud<pcl::PointXYZ>::Ptr ann_cloud( new pcl::PointCloud<pcl::PointXYZ>() );
    {
        ann_cloud->reserve( points.size() );
        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            pcl::PointXYZ pnt;
            pnt.x = points[pid].template pos()(0); // convert PointPrimitive to Eigen::Matrix, and get (0)
            pnt.y = points[pid].template pos()(1);
            pnt.z = points[pid].template pos()(2);
            ann_cloud->push_back( pnt );
        }
    }
    typename pcl::search::KdTree<pcl::PointXYZ>::Ptr tree( new pcl::search::KdTree<pcl::PointXYZ> );
    tree->setInputCloud( ann_cloud );

    Patches patches; patches.reserve( std::max(1.5*sqrt(points.size()),10.) );

    // get unassigned point
    std::vector<bool> assigned( points.size(), false );
    std::vector<bool> visited( points.size(), false );

    std::vector<float>  sqr_dists( nn_K );
    std::vector< int >  neighs( nn_K );
    int                 found_points_count  = 0;
    pcl::PointXYZ       searchPoint;
    const _Scalar       max_dist            = patchPatchDistanceFunctor.getSpatialThreshold() * _Scalar(3.5); // longest axis of ellipse)

    // look for neighbours, merge most similar
    while ( starting_cands.size() )
    {
        // remove point from unassigned
        int pid = starting_cands.front();
        starting_cands.pop_front();

        if ( visited[pid] && !assigned[pid] )
            std::cerr << "[" << __func__ << "]: " << " visited but not assigned...:JSDL:FJSD:LKFSDK:SFDJ" << std::endl;

        if ( visited[pid] )
            continue;

        // add to new cluster
        if ( !assigned[pid] )
        {
            PatchT tmp_patch; tmp_patch.push_back( segmentation::PidLid(pid,-1) );
            patches.push_back( tmp_patch );
            patches.back().update( points );
            assigned[ pid ] = true;
        }

        // look for unassigned neighbours
        searchPoint.x = points[ pid ].template pos()(0);
        searchPoint.y = points[ pid ].template pos()(1);
        searchPoint.z = points[ pid ].template pos()(2);
        found_points_count = tree->radiusSearch( searchPoint, max_dist, neighs, sqr_dists, 0);

        if ( neighs.size() != found_points_count )
            std::cerr << "as;dlkjfa;lkdsjfl;asdjkf neighs.size() " << neighs.size() << " != " << found_points_count << " found_points_count" << std::endl;

        for ( size_t pid_id = 0; pid_id != neighs.size(); ++pid_id )
        {
            const int pid2 = neighs[ pid_id ];
            if ( !assigned[pid2] )
            {
                //_Scalar diff = GF2::angleInRad( patches.back().template dir(), points[pid2].template dir() );

                // location from point, but direction is the representative's
                _PointPrimitiveT p1_proxy(points[pid].template pos(), patches.back().template dir());
                PatchT p2_proxy; p2_proxy.push_back( segmentation::PidLid(pid2,-1) ); p2_proxy.update( points );

                //if ( (sqrt(sqr_dists[pid_id]) < spatial_thresh) && (diff < ang_thresh) )                                      // original condition
                if ( patchPatchDistanceFunctor.template eval<_PointPrimitiveT>(p1_proxy, p2_proxy, points, NULL) < patchPatchDistanceFunctor.getThreshold() )
                {
                    patches.back().push_back( segmentation::PidLid(pid2,-1) );
                    patches.back().updateWithPoint( points[pid2] );
                    assigned[pid2] = true;
                    // enqueue for visit
                    starting_cands.push_front( pid2 );
                }
            }
        }

        visited[pid] = true;
    }

    // copy patches to groups
    groups_arg.insert( groups_arg.end(), patches.begin(), patches.end() );

    // assign points to patches
    _tagPointsFromGroups<_PointPrimitiveT,_Scalar>
                        ( points, groups_arg, patchPatchDistanceFunctor, gid_tag_name );

    return EXIT_SUCCESS;
} // ...Segmentation::regionGrow()

/*  \brief                  Step 1. Generates primitives from a cloud. Reads "cloud.ply" and saves "candidates.txt".
 *  \param argc             Contains --cloud cloud.ply, and --scale scale.
 *  \param argv             Contains --cloud cloud.ply, and --scale scale.
 *  \return                 EXIT_SUCCESS.
 *  \post                   "patches.txt" and "points_primitives.txt" on disk in "cloud.ply"'s parent path.
 */
template < class _PrimitiveT
         , class _PrimitiveContainerT
         , class _PointPrimitiveT
         , class _PointContainerT
         , typename _Scalar
         >
int
Segmentation::segmentCli( int    argc
                        , char** argv )
{
    int err = EXIT_SUCCESS;

    CandidateGeneratorParams<_Scalar> generatorParams;
    std::string                 cloud_path              = "./cloud.ply";
    _Scalar                     angle_gen               = M_PI_2;
    std::string                 mode_string             = "representative_sqr";
    std::vector<std::string>    mode_opts               = { "representative_sqr" };

    // parse input
    if ( err == EXIT_SUCCESS )
    {
        bool valid_input = true;

        // scale
        if ( (pcl::console::parse_argument( argc, argv, "--scale", generatorParams.scale) < 0) )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale is compulsory" << std::endl;
            valid_input = false;
        }

        // cloud
        if ( (pcl::console::parse_argument( argc, argv, "--cloud", cloud_path) < 0)
             && !boost::filesystem::exists( cloud_path ) )
        {
            std::cerr << "[" << __func__ << "]: " << "--cloud does not exist: " << cloud_path << std::endl;
            valid_input = false;
        }

        pcl::console::parse_argument( argc, argv, "--angle-limit", generatorParams.angle_limit );
        pcl::console::parse_argument( argc, argv, "--patch-dist-limit", generatorParams.patch_dist_limit_mult ); // gets multiplied by scale
        pcl::console::parse_argument( argc, argv, "--mode", mode_string );
        generatorParams.parsePatchDistMode( mode_string );

        if ( pcl::console::find_switch( argc, argv, "--patch-refit" ) )
        {
            std::cerr << "[" << __func__ << "]: " << "--patch-refit option has been DEPRECATED. exiting." << std::endl;
            return EXIT_FAILURE;
        }
        pcl::console::parse_argument( argc, argv, "--angle-gen", angle_gen );

        // print usage
        {
            std::cerr << "[" << __func__ << "]: " << "Usage:\t " << argv[0] << " --segment \n";
            std::cerr << "\t --cloud " << cloud_path << "\n";
            std::cerr << "\t --scale " << generatorParams.scale      << "\n";

            // linkage mode (full_min, full_max, squared_min, repr_min)
            std::cerr << "\t [--mode *" << generatorParams.printPatchDistMode() << "*\t";
            for ( size_t m = 0; m != mode_opts.size(); ++m )
                std::cerr << "|" << mode_opts[m];
            std::cerr << "]\n";

            std::cerr << "\t [--angle-limit " << generatorParams.angle_limit << "]\n";
            std::cerr << "\t [--patch-dist-limit " << generatorParams.patch_dist_limit_mult << "]\n";
            std::cerr << "\t [--angle-gen " << angle_gen << "]\n";
            std::cerr << std::endl;

            if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
                return EXIT_FAILURE;
        }

        if ( boost::filesystem::is_directory(cloud_path) )
        {
            cloud_path += "/cloud.ply";
        }

        if ( !boost::filesystem::exists(cloud_path) )
        {
            std::cerr << "[" << __func__ << "]: " << "cloud file does not exist! " << cloud_path << std::endl;
            return EXIT_FAILURE;
        }
    } // ... parse input

    // Read desired angles
    if ( EXIT_SUCCESS == err )
    {
        processing::appendAnglesFromGenerator( generatorParams.angles, angle_gen, true );
    } //...read angles

    // Read points
    _PointContainerT points;
    if ( EXIT_SUCCESS == err )
    {
        err = io::readPoints<_PointPrimitiveT>( points, cloud_path );
        if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
    } //...read points

    //_____________________WORK_______________________
    //_______________________________________________

    // orientPoints
    if ( EXIT_SUCCESS == err )
    {
        err = Segmentation::orientPoints<_PointPrimitiveT,_PrimitiveT>( points, generatorParams.scale, generatorParams.nn_K );
        if ( err != EXIT_SUCCESS ) std::cerr << "[" << __func__ << "]: " << "orientPoints exited with error! Code: " << err << std::endl;
    } //...orientPoints

    _PrimitiveContainerT initial_primitives;
    if ( EXIT_SUCCESS == err )
    {
        switch ( generatorParams.patch_dist_mode )
        {
            case CandidateGeneratorParams<_Scalar>::REPR_SQR:
            {
                // "representative min:" merge closest representative angles, IF smallest spatial distance between points < scale * patch_dist_limit.
                RepresentativeSqrPatchPatchDistanceFunctorT< _Scalar,SpatialPatchPatchSingleDistanceFunctorT<_Scalar>
                                                        > patchPatchDistanceFunctor( generatorParams.scale * generatorParams.patch_dist_limit_mult
                                                                                   , generatorParams.angle_limit
                                                                                   , generatorParams.scale
                                                                                   , generatorParams.patch_spatial_weight );
                err = Segmentation::patchify<_PrimitiveT>( initial_primitives    // tagged lines at GID with patch_id
                                            , points                            // filled points with directions and tagged at GID with patch_id
                                            , generatorParams.scale
                                            , generatorParams.angles
                                            , patchPatchDistanceFunctor
                                            , generatorParams.nn_K
                                            );
            }
                break;

            default:
                std::cerr << "unknown patch patch distance mode!" << std::endl;
                err = EXIT_FAILURE;
                break;
        }

        if ( err != EXIT_SUCCESS ) std::cerr << "[" << __func__ << "]: " << "patchify exited with error! Code: " << err << std::endl;
    }

    // Save point GID tags
    if ( EXIT_SUCCESS == err )
    {
        std::string assoc_path = boost::filesystem::path( cloud_path ).parent_path().string() + "/" + "points_primitives.csv";

        util::saveBackup( assoc_path );
        err = io::writeAssociations<_PointPrimitiveT>( points, assoc_path );

        if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "saveBackup or writeAssociations exited with error! Code: " << err << std::endl;
        else                        std::cout << "[" << __func__ << "]: " << "wrote to " << assoc_path << std::endl;

    } //...save Associations

    // save primitives
    if ( EXIT_SUCCESS == err )
    {
        std::string candidates_path = boost::filesystem::path( cloud_path ).parent_path().string() + "/" + "patches.csv";

        util::saveBackup( candidates_path );
        err = io::savePrimitives<_PrimitiveT,typename _PrimitiveContainerT::value_type::const_iterator>( /* what: */ initial_primitives, /* where_to: */ candidates_path );

        if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "saveBackup or savePrimitive exited with error! Code: " << err << std::endl;
        else                        std::cout << "[" << __func__ << "]: " << "wrote to " << candidates_path << std::endl;
    } //...save primitives

    return err;
} // ...Segmentation::segmentCli()

template < class    _PointPrimitiveT
         , typename _Scalar
         , class    _PointPatchDistanceFunctorT
         , class    _PatchT
         , class    _PointContainerT
         >  int
Segmentation::_tagPointsFromGroups( _PointContainerT                 & points
                                  , _PatchT                     const& groups
                                  , _PointPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                  , int                         const  gid_tag_name )
{
    // set group tag for grouped points
    std::vector<bool> visited( points.size(), false );
    for ( int gid = 0; gid != groups.size(); ++gid )
        for ( size_t pid_id = 0; pid_id != groups[gid].size(); ++pid_id )
        {
            const int pid = groups[gid][pid_id].first;

            visited.at(pid) = true; // TODO: change to []
            points[pid].setTag( gid_tag_name, gid );
        }

    // add left out points to closest patch
    for ( int pid = 0; pid != points.size(); ++pid )
    {
        if ( visited[pid] )
            continue;

        _Scalar     min_dist = std::numeric_limits<_Scalar>::max(), dist;
        int         min_gid  = 0;
        for ( int gid = 1; gid != groups.size(); ++gid )
        {
            if ( (dist = pointPatchDistanceFunctor.template evalSpatial<_PointPrimitiveT>(pid, groups[gid], points)) < min_dist )
            {
                min_dist = dist;
                min_gid  = gid;
            }
        }

        points[pid].setTag( gid_tag_name, min_gid );
    }

    return EXIT_SUCCESS;
} // ...Segmentation::tagPointsFromGroups()

} //...ns GF2

#endif // GF2_SEGMENTATION_HPP
