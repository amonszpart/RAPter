#ifndef GF2_MERGING_HPP
#define GF2_MERGING_HPP

#include "globfit2/optimization/merging.h"

#include "globfit2/parameters.h"
#include "globfit2/visualization/visualization.h"
#include "globfit2/io/io.h"
#include "globfit2/processing/util.hpp"          //getPopulations()
#include "globfit2/optimization/patchDistanceFunctors.h" // RepresentativeSqrPatchPatchDistanceFunctorT
#include "globfit2/util/util.hpp"

#define CHECK(err,text) { if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << text << " returned an error! Code: " << err << std::endl; }

namespace GF2 {

template < class    _PrimitiveContainerT
         , class    _PointContainerT
         , typename _Scalar
         , class    _PointPrimitiveT
         , class    _PrimitiveT
         >
int
Merging::mergeCli( int argc, char** argv )
{
    MergeParams<_Scalar> params;

    std::string cloud_path = "cloud.ply",
                prims_path = "primitives.bonmin.csv",
                assoc_path = "points_primitives.csv";
    _Scalar     angle_gen  = M_PI_2;
    // parse params
    {
        bool valid_input = true;

        valid_input &= pcl::console::parse_argument( argc, argv, "--scale", params.scale ) >= 0;
        valid_input &= pcl::console::parse_argument( argc, argv, "--prims", prims_path   ) >= 0;

        // cloud
        pcl::console::parse_argument( argc, argv, "--cloud", cloud_path );
        valid_input &= boost::filesystem::exists( cloud_path );

        pcl::console::parse_argument( argc, argv, "--angle-gen", angle_gen );
        pcl::console::parse_argument( argc, argv, "--adopt", params.do_adopt );
        pcl::console::parse_argument( argc, argv, "--thresh-mult", params.spatial_threshold_mult );
        pcl::console::parse_argument( argc, argv, "--assoc", assoc_path );
        pcl::console::parse_argument( argc, argv, "-a", assoc_path );

        std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt --formulate\n"
                  << "\t--scale " << params.scale << "\n"
                  << "\t--prims " << prims_path << "\n"
                  << "\t--cloud " << cloud_path << "\n"
                  << "\t-a,--assoc " << assoc_path << "\n"
                  << "\t[--angle-gen " << angle_gen << "]\n"
                  << "\t[--adopt " << params.do_adopt << "]\n"
                  << "\t[--thresh-mult " << params.spatial_threshold_mult << "]\n"
                  << std::endl;

        if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale, --prims are compulsory, --cloud needs to exist" << std::endl;
            return EXIT_FAILURE;
        }
    } // ... parse params

    // read primitives
    typedef std::vector<_PrimitiveT>    PatchT;
    typedef std::map   < int, PatchT >  PrimitiveMapT;
    _PrimitiveContainerT prims;
    PrimitiveMapT        prims_map;
    {
        io::readPrimitives<_PrimitiveT, PatchT>( prims, prims_path, &prims_map );
    } //...read primitives

    // Read points
    _PointContainerT points;
    {
        io::readPoints<_PointPrimitiveT>( points, cloud_path );
    }

    // Read desired angles
    processing::appendAnglesFromGenerator( params.angles, angle_gen, true );

    // associations
    std::vector<std::pair<int,int> > points_primitives;
    io::readAssociations( points_primitives, assoc_path, NULL );
    for ( size_t i = 0; i != points.size(); ++i )
    {
        // error check
        if ( i > points_primitives.size() )
        {
            std::cerr << "more points than associations..." << std::endl;
            return EXIT_FAILURE;
        }

        // store association in point
        points[i].setTag( _PointPrimitiveT::GID, points_primitives[i].first );

        // error check 2
        if ( points[i].getTag(_PointPrimitiveT::GID) == -1 )
            std::cerr << "[" << __func__ << "]: " << "point assigned to patch with id -1" << std::endl;
    }

    //____________________________WORK____________________________

    // ADOPT
    if ( params.do_adopt )
    {
        adoptPoints<GF2::MyPointPrimitiveDistanceFunctor, _PointPrimitiveT, _PrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>
                    ( points, prims_map, params.scale, params.do_adopt );
    }

    // MERGE
    std::cout << "starting mergeSameDirGids" << std::endl; fflush(stdout);
    RepresentativeSqrPatchPatchDistanceFunctorT<_Scalar, SpatialPatchPatchSingleDistanceFunctorT<_Scalar> >
            patchPatchDistFunct( params.scale * params.patch_dist_limit_mult
                               , params.angle_limit
                               , params.scale
                               , params.patch_spatial_weight );
    PrimitiveMapT out_prims;
    mergeSameDirGids<_PrimitiveT, _PointPrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>
            ( out_prims, points, prims_map, params.scale, params.spatial_threshold_mult * params.scale, params.parallel_limit, patchPatchDistFunct );

    // SAVE
    std::string o_path;
    int         iteration = 0;
    {
        iteration = util::parseIteration( prims_path );
        std::stringstream ss;
        size_t it_loc = prims_path.find("_it");
        std::string fname = prims_path.substr( 0, it_loc );
        ss << fname << "_merged_it" << iteration << ".csv";
        o_path = ss.str();
        io::savePrimitives   <_PrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>( out_prims, o_path );
        std::cout << "wrote " << o_path << std::endl;
    }

    {
        std::stringstream ss;
        ss << "points_primitives_it" << iteration << ".csv";
        io::writeAssociations<_PointPrimitiveT>( points, ss.str() );
        std::cout << "wrote " << ss.str() << std::endl;
    }

    std::cout << "stopped mergeSameDirGids" << std::endl; fflush(stdout);
    return EXIT_SUCCESS;
}//...Merging::mergeCli()

/*! \brief Greedily assigns points with GID-s that are not in prims to prims that explain them.
*        Unambiguous assignments go through first, than based on proximity, capped by scale.
*
* \tparam _PointPrimitiveDistanceFunctor Has an eval function for a point and all primitives, to calculate the distance from point to primitive. Concept: \ref MyPointPrimitiveDistanceFunctor.
* \tparam _PointPrimitiveT     Wraps a point, exposing pos() and dir() functions. Concept: \ref GF2::PointPrimitive.
* \tparam _PrimitiveT          Wraps a primitive, exposing pos() and dir() functions. Concept: \ref GF2::LinePrimitive2.
* \tparam _PointContainerT     Holds the points. Concept: std::vector< \ref GF2::PointPrimitive >.
* \tparam _PrimitiveContainerT Holds the primitives grouped by GID in a two level structure. Concept: std::map< int, std::vector< \ref GF2::LinePrimitive2 > >.
* \tparam _Scalar              Floating point precision of primitives, points, etc. Concept: \ref GF2::PointPrimitive::Scalar.
* \param[in,out] points        Contains the points, some assigned, some to be assigned to the primitives in prims.
* \param[in]     prims         Contains some primitives tagged with GID and DIR_GID. GID defines the assignment between points and primitives.
* \param[in]     scale         Distance threshold parameter.
* \param[in]     mode          1: re-assign un-ambiguous points (1 adopter); 2: first re-assign unambiguous, then closest, if inside explaining primitive's scale.
*/
template < class _PointPrimitiveDistanceFunctor
         , class _PointPrimitiveT
         , class _PrimitiveT
         , class _inner_const_iterator
         , class _PointContainerT
         , class _PrimitiveContainerT
         , typename _Scalar >
int Merging::adoptPoints( _PointContainerT          & points
                        , _PrimitiveContainerT const& prims
                        , _Scalar              const  scale
                        , char                 const  mode )
{
    //typedef GF2::MyPointPrimitiveDistanceFunctor _PointPrimitiveDistanceFunctor;
    typedef typename _PrimitiveContainerT::const_iterator             outer_const_iterator;
    //typedef typename outer_const_iterator::value_type::const_iterator inner_const_iterator;

    int err = EXIT_SUCCESS;

    // select unassigned points
    std::deque<int> orphan_pids;
    {
        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            int gid = points[pid].getTag( _PointPrimitiveT::GID );
            if (    (prims.find(gid) == prims.end())
                 || (!containers::valueOf<_PrimitiveT>(prims.find(gid)).size()) )
            //if ( !prims.at(gid).size() ) // no primitives with this GID
            {
                points[pid].setTag   ( _PointPrimitiveT::GID, -2 );
                orphan_pids.push_back( pid );
            }
        }
    }

    int change = 0, iteration = 0;
    do
    {
        std::cout << "[" << __func__ << "]: " << "iteration " << iteration << std::endl;
        // reset re-assignment counter
        change = 0;

        // store possible primitives to assign points[pid] to
        typedef std::pair<int,int> GidLid1;
        std::map<int, std::vector<GidLid1> > adopter_gids; // < pid, [ <gid,lid1>, ... ] >
        {
            // iterate over unassigned points
            std::deque<int>::const_iterator points_it_end = orphan_pids.end();
            for ( std::deque<int>::const_iterator points_it = orphan_pids.begin(); points_it != points_it_end; ++points_it )
            {
                // store point id
                const int pid = *points_it;

                // for patches
                //typename PrimitiveMapT::const_iterator end_it = prims.end();
                outer_const_iterator end_it = prims.end();
                for ( outer_const_iterator it = prims.begin(); it != end_it; ++it )
                {
                    // record primitive id in patch
                    int lid1 = 0;
                    // for primitives in patch
                    _inner_const_iterator end_it2 = it->second.end();
                    for ( _inner_const_iterator it2 = it->second.begin(); it2 != end_it2; ++it2, ++lid1 )
                    {
                        // distance from point to primitive
                        _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], *it2 );
                        // store if inside scale
                        if ( dist < scale )
                        {
                            adopter_gids[pid].push_back( std::pair<int,int>(it->first,lid1) );
                        } //...if dist

                    } //...for prims_in_patch
                } //...for patches
            } //...for points
        } //...adopters

        // reassign points
        if ( iteration == 0 ) //
        {
            for ( std::deque<int>::iterator points_it = orphan_pids.begin(); points_it != orphan_pids.end(); /*nothing*/ )
            {
                // store point id
                const int pid = *points_it;
                //std::cout << "[" << __func__ << "]: " << "point " << pid << " have " << adopter_gids[pid].size() << " adopters" << std::endl;

                if ( !adopter_gids[pid].size() ) std::cerr << "[" << __func__ << "]: " << "point " << pid << " seems to be an outlier...no points could adopt it..." << std::endl;

                // if point has only one possible reassignment - do assign it
                if ( adopter_gids[pid].size() == 1 )
                {
                    points[pid].setTag( _PointPrimitiveT::GID, adopter_gids[pid][0].first );
                    points_it = orphan_pids.erase( points_it ); // returns iterator to next element
                    ++change;
                }
                else
                {
                    ++points_it; // increment only, if no erase has been done
                }
            }
        }
        else if ( mode == 2 )
        {
            int closest_gid = -3; // gid of point, that is closest to an orphan, termination crit
            do
            {
                // reset
                closest_gid = -3;
                _Scalar                   closest_distance  = std::numeric_limits<_Scalar>::max();
                std::deque<int>::iterator closest_orphan    = orphan_pids.begin();

                for ( std::deque<int>::iterator points_it = orphan_pids.begin(); points_it != orphan_pids.end(); ++points_it )
                {
                    // store point id
                    const int pid = *points_it; // orphan_pid

                    // look for the closest assigned point from possible adopters
                    for ( int pid2 = 0; pid2 != points.size(); ++pid2 )
                    {
                        // rule out itself
                        if ( pid2 == pid ) continue;
                        // rule out other orphans
                        const int gid2 = points[pid2].getTag(_PointPrimitiveT::GID);
                        if ( gid2 < 0 )  continue;

                        _Scalar dist = (points[pid].template pos() - points[pid2].template pos()).norm();
                        // store if closer, than earlier orphan, and if inside scale to given primitive (can be explained)
                        if ( (dist < closest_distance) && _PointPrimitiveDistanceFunctor::template eval<_Scalar>(points[pid], prims.at(gid2).at(0)) )
                        {
                            if ( prims.at(gid2).size() > 1 ) std::cerr << "[" << __func__ << "]: " << "two lines in one patch, not prepared here, fix this" << std::endl;

                            closest_distance = dist;
                            closest_gid      = points[pid2].getTag( _PointPrimitiveT::GID );
                            closest_orphan   = points_it;
                        } //...if dist <
                    } //...for cloud
                } //...for all orphans

                // check, if any luck
                if ( closest_gid >= 0 )
                {
                    points[ *closest_orphan ].setTag( _PointPrimitiveT::GID, closest_gid );
                    std::cout << "erasing orphan " << std::distance( orphan_pids.begin(), closest_orphan ) << "/" << orphan_pids.size();
                    orphan_pids.erase( closest_orphan );
                    std::cout << " size now " << orphan_pids.size() << std::endl;
                    ++change;
                }
            } while ( (closest_gid >= 0) && orphan_pids.size() );

        } //...if iteration != 0

        ++iteration;
    } while ( change ); // stop, if no new points were reassigned

    return err;
} //...adoptPoints()

/*! \brief Decides, if two patches are adjacent and have the same direction. Used in \ref erging::mergeSameDirGids().
 *
 */
//template <typename _Scalar>
//inline bool decide_merge( _Scalar min_dist, _Scalar threshold, _Scalar angle, _Scalar parallel_limit )
//{
//    return (min_dist < threshold) && (angle < parallel_limit);
//} //...decide_merge

template <class _PrimitiveT, typename _Scalar>
//inline _Scalar segmentDistance(
inline bool decide_merge(
        std::vector<Eigen::Matrix<_Scalar,3,1> > const& extrema0
        , _PrimitiveT const& l0
        , std::vector<Eigen::Matrix<_Scalar,3,1> > const& extrema1
        , _PrimitiveT const& l1
        , _Scalar scale)
{

    // We have two lines l0 and l1, respectively defined by
    // l0a-l0b and l1a-l1b
    //  1. We project l0a and l0b to l1 (orthogonally to l1), and compute the norm of l0a-proj(l0a,l1) and l0b-proj(l0b,l1)
    //  2. Compute the other way
    //  2. The two lines can be merged both norms of one of the two ways are below the scale parameters
    const Eigen::Matrix<_Scalar,3,1> & l0a = extrema0[0];
    const Eigen::Matrix<_Scalar,3,1> & l0b = extrema0[1];
    const Eigen::Matrix<_Scalar,3,1> & l1a = extrema1[0];
    const Eigen::Matrix<_Scalar,3,1> & l1b = extrema1[1];

    Eigen::Matrix<_Scalar,3,1> n0 = l0.template normal<_Scalar>();
    Eigen::Matrix<_Scalar,3,1> n1 = l1.template normal<_Scalar>();


    return ( std::abs(n1.dot(l0a-l1a)) <= scale && // get l0a-proj(l0a,l1)
             std::abs(n1.dot(l0b-l1a)) <= scale ) ||
           ( std::abs(n0.dot(l1a-l0a)) <= scale && // get l0a-proj(l0a,l1)
             std::abs(n0.dot(l1b-l0a)) <= scale )
            ;


    /*_Scalar min_dist = std::numeric_limits<_Scalar>::max();

    for ( int i = 0; i != extrema0.size(); ++i )
        for ( int j = 0; j != extrema1.size(); ++j )
        {
            Eigen::Matrix<_Scalar,3,1> const& e0 = extrema0[i];
            Eigen::Matrix<_Scalar,3,1> const& e1 = extrema1[j];

            Eigen::Matrix<_Scalar,3,1> p0 = l1.projectPoint( e0 );
            Eigen::Matrix<_Scalar,3,1> p1 = l0.projectPoint( e1 );

            _Scalar paramCoord0 = (p0 - extrema1[0]).norm() / (extrema1[1]-extrema1[0]).norm();
            _Scalar paramCoord1 = (p1 - extrema0[0]).norm() / (extrema0[1]-extrema0[0]).norm();

            std::cout << "paramCoord0: " << paramCoord0
                      << ", paramCoord1: " << paramCoord1 << std::endl;

            bool in0 = (paramCoord0 >= _Scalar(0.) && paramCoord0 <= _Scalar(1.) );
            bool in1 = (paramCoord1 >= _Scalar(0.) && paramCoord1 <= _Scalar(1.) );
            switch ( in0 + in1 )
            {
                case 0: std::cout  << "non";
                    min_dist = std::min( min_dist, (e1 - e0).norm() );
                    break;
                case 1: std::cout  << "one";
                         if ( in0 ) { min_dist = std::min( min_dist, (p1-e1).norm()); std::cout << " (p1-e1).norm(): " << (p1-e1).norm(); }
                    else if ( in1 ) { min_dist = std::min( min_dist, (p0-e0).norm()); std::cout << " (p0-e0).norm(): " << (p0-e0).norm(); }
                    break;
                case 2:
                    std::cout  << "mutual std::min( " << min_dist << ", std::min( " << (p1-e1).norm() << ", " << (p0-e0).norm() << " ))" << std::endl;
                    min_dist = std::min( min_dist, std::min( (p1-e1).norm(), (p0-e0).norm()) );
                    break;
            }
            std::cout << std::endl;
        }

    return min_dist;*/
}

/*! \brief Merges adjacent patches that have the same direction ID or are almost parallel.
 *
 *  \tparam _PatchPatchDistanceFunctorT  Concept: \ref GF2::RepresentativeSqrPatchPatchDistanceFunctorT.
 *  \param[in] patchPatchDistFunct       Distance functor between two patches, to define adjacency.
 *  \param[in] spatial_threshold         Two extrema should be at least this close to be merged. Concept: \ref MergeParams::spatial_threshold_mult == 3 * scale.
 *  \sa #decide_merge()
 */
template < class    _PrimitiveT
         , class    _PointPrimitiveT
         , class    _inner_const_iterator
         , class    _PrimitiveContainerT
         , class    _PointContainerT
         , typename _Scalar
         , class    _PatchPatchDistanceFunctorT>
int Merging::mergeSameDirGids( _PrimitiveContainerT             & out_primitives
                             , _PointContainerT                 & points
                             , _PrimitiveContainerT        const& primitives
                             , _Scalar                     const  scale
                             , _Scalar                     const  spatial_threshold
                             , _Scalar                     const  parallel_limit
                             , _PatchPatchDistanceFunctorT const& patchPatchDistFunct )
{
    typedef typename _PrimitiveContainerT::const_iterator      outer_const_iterator;
    typedef           std::vector<Eigen::Matrix<_Scalar,3,1> > ExtremaT;
    typedef           std::map   < int, ExtremaT>              LidExtremaT;
    typedef           std::map   < int, LidExtremaT >          GidLidExtremaT;
    typedef           std::pair  < int, int> GidLid;

    int err = EXIT_SUCCESS;

    // Populations
    GidPidVectorMap populations; // populations[gid] == std::vector<int> {pid0,pid1,...}
    if ( EXIT_SUCCESS == err )
    {
        err = processing::getPopulations( populations, points );
        CHECK( err, "getPopulations" );
    }

    // Extrema
    GidLidExtremaT extrema; // <gid,lid> -> vector<x0, x1, ...>
    if ( EXIT_SUCCESS == err )
    {
        // for all patches
        for ( outer_const_iterator outer_it  = primitives.begin();
                                  (outer_it != primitives.end())  && (EXIT_SUCCESS == err);
                                 ++outer_it )
        {
            int gid  = -2; // (-1 is "unset")
            int lid = 0; // linear index of primitive in container (to keep track)
            // for all directions
            for ( _inner_const_iterator inner_it  = containers::valueOf<_PrimitiveT>(outer_it).begin();
                                       (inner_it != containers::valueOf<_PrimitiveT>(outer_it).end()) && (EXIT_SUCCESS == err);
                                      ++inner_it, ++lid )
            {
                // save patch gid
                if ( gid == -2 )
                {
                    gid = inner_it->getTag( _PrimitiveT::GID );
                    // sanity check
                    if ( extrema.find(gid) != extrema.end() )
                    {
                        std::cerr << "[" << __func__ << "]: " << "GID not unique for patch...:-S" << std::endl;
                    }
                }

                err = inner_it->template getExtent<_PointPrimitiveT>
                                                           ( extrema[gid][lid]
                                                           , &points
                                                           , scale
                                                           , populations[gid].size() ? &(populations[gid]) : NULL );
            } //...for primitives
        } //...for patches

        CHECK( err, "calcExtrema" )
    } //...getExtrema

    std::cout << "SCALE = " << scale << std::endl;

    typedef std::map< GidLid, GidLid > AliasesT;
    AliasesT aliases;
    // outer traversal
    typename GidLidExtremaT::const_iterator gid_end_it = extrema.end();
    for ( typename GidLidExtremaT::const_iterator gid_it = extrema.begin(); gid_it != gid_end_it; ++gid_it )
    {
        typename LidExtremaT::const_iterator prim_end_it = gid_it->second.end();
        for ( typename LidExtremaT::const_iterator prim_it = gid_it->second.begin(); prim_it != prim_end_it; ++prim_it )
        {
            // inner traversal
            typename GidLidExtremaT::const_iterator gid_end_it1 = extrema.end();
            for ( typename GidLidExtremaT::const_iterator gid_it1 = extrema.begin(); gid_it1 != gid_end_it1; ++gid_it1 )
            {
                // offset to start from same as outer (don't want to go beyond, so not +1, but "continue" later)
                if ( gid_it1 == extrema.begin() )
                    std::advance( gid_it1, std::distance<typename GidLidExtremaT::const_iterator>(extrema.begin(),gid_it) );

                typename LidExtremaT::const_iterator prim_end_it1 = gid_it1->second.end();
                for ( typename LidExtremaT::const_iterator prim_it1 = gid_it1->second.begin(); prim_it1 != prim_end_it1; ++prim_it1 )
                {
                    // offset to start from same as outer (don't want to go beyond, so not +1, but "continue" if exactly the same)
                    if ( prim_it1 == gid_it1->second.begin() )
                        std::advance( prim_it1, std::distance<typename LidExtremaT::const_iterator>(gid_it->second.begin(),prim_it) );

                    // skip, if itself
                    if ( (gid_it == gid_it1) && (prim_it == prim_it1) )
                        continue;

                    int gid0 = gid_it->first ,
                        gid1 = gid_it1->first;
                    int lid0 = std::distance<typename LidExtremaT::const_iterator>( gid_it ->second.begin(), prim_it ),
                        lid1 = std::distance<typename LidExtremaT::const_iterator>( gid_it1->second.begin(), prim_it1 );

                    // get minimum endpoint distance
                    ExtremaT const& extrema0 = prim_it->second,
                                    extrema1 = prim_it1->second;
                    _Scalar min_dist = std::numeric_limits<_Scalar>::max();
#if 1
                    _PrimitiveT const& prim0 = primitives.at(gid0).at(lid0);
                    _PrimitiveT const& prim1 = primitives.at(gid1).at(lid1);
                    //min_dist = segmentDistance( extrema0, prim0, extrema1, prim1 );
#else
                    for ( int i = 0; i != extrema0.size(); ++i )
                        for ( int j = 0; j != extrema1.size(); ++j )
                        {
                            _Scalar dist = (extrema0[i] - extrema1[j]).norm();
                            if ( dist < min_dist )
                            {
                                min_dist = dist;
                            }
                        }
#endif


                    int dir0 = primitives.at( gid0 ).at( lid0 ).getTag( _PrimitiveT::DIR_GID ); //used for print only
                    int dir1 = primitives.at( gid1 ).at( lid1 ).getTag( _PrimitiveT::DIR_GID ); //used for print only
                    //_Scalar ang = angleInRad( primitives.at( gid0 ).at( lid0 ).template dir(),  //used for print only
                    //                          primitives.at( gid1 ).at( lid1 ).template dir() );

                    std::cout << "testing (" << gid0 << "," << dir0 << ")" << " vs. " << "(" << gid1 << "," << dir1 << ")\t" << endl;
                    //          << min_dist << " < " <<  spatial_threshold << " && " << ang << " < " << parallel_limit;
                    //if ( decide_merge(min_dist, spatial_threshold, ang, parallel_limit) ) // (min_dist < 3. * scale) && (dir0 == dir1) && (ang < parallel_limit) )
                    if (decide_merge( extrema0, prim0, extrema1, prim1, scale ))
                    {
                        std::cout << "would merge "
                                  << "(" << gid0 << "," << lid0 << "," << dir0 << ")"
                                  << " with "
                                  << "(" << gid1 << "," << lid1 << "," << dir1 << ")"
                                  << std::endl;
                        fflush( stdout );

                        GidLid key0( gid0, lid0 ),
                               key1( gid1, lid1 );

                        AliasesT::const_iterator alias1 = aliases.find( key1 ), alias0;
                        if ( alias1 != aliases.end() )
                            aliases[ key0 ] = alias1->second; // if key1 is already merged with *key1, then key0 should be merged with *key1 as well.
                        else if ( (alias0 = aliases.find(key0)) != aliases.end() )
                            aliases[ key1 ] = alias0->second;  // if key0 is already merged with *key0, then key1 should be merged with *key0 as well.
                        else
                            aliases[ key1 ] = key0;
                    }
                    else
                        std::cout << " NO" << std::endl;
                }
            }
        }
    }

    for ( auto it = aliases.begin(); it != aliases.end(); ++it )
    {
        std::cout << "(" << it->first.first << "," << it->first.second << ") - (" << it->second.first << "," << it->second.second << ")" << std::endl;
    }


    for ( outer_const_iterator outer_it  = primitives.begin();
                              (outer_it != primitives.end())  && (EXIT_SUCCESS == err);
                             ++outer_it )
    {
        int gid = -2; // (-1 is "unset")
        int lid =  0; // linear index of primitive in container (to keep track)
        // for all directions
        for ( _inner_const_iterator inner_it  = containers::valueOf<_PrimitiveT>(outer_it).begin();
                                   (inner_it != containers::valueOf<_PrimitiveT>(outer_it).end()) && (EXIT_SUCCESS == err);
                                  ++inner_it, ++lid )
        {
            if ( gid == -2 )
                gid = inner_it->getTag( _PrimitiveT::GID );
            GidLid key( gid, lid );
            if ( aliases.find(key) == aliases.end() )
            {
                containers::add<_PrimitiveT>( out_primitives, gid, *inner_it );
#warning TODO: average instead of throw away
            }
        }
    }

    for ( int pid = 0; pid != points.size(); ++pid )
    {
        const int gid = points[pid].getTag( _PointPrimitiveT::GID );

        // find an entry with an alias for this gid
        AliasesT::const_iterator alias_it = std::find_if( aliases.begin(), aliases.end(), [&gid]( AliasesT::value_type const& pair ) { return pair.first.first == gid; } );
        if ( alias_it != aliases.end() )
        {
            std::cout << "setting gid " << gid << " to " << alias_it->second.first << std::endl;
            points[pid].setTag( _PointPrimitiveT::GID, alias_it->second.first );
        }
    }

    //patchPatchDistanceFunctor.template eval<_PointPrimitiveT>(p1_proxy, p2_proxy, points, NULL)

    // debug
//    std::ofstream dbg( "lines.plot" );
//    for ( auto it = extrema.begin(); it != extrema.end(); ++it )
//    {
//        for ( int i = 0; i != it->second.size(); ++i )
//        {
//            dbg << extrema[it->first][i][0].transpose() << "\n"
//                << extrema[it->first][i][1].transpose() << "\n\n";
//        }
//    }
//    dbg.close();
    return err;
} //...Merging::mergeSameDirGids()

} //...namespace GF2

#undef CHECK

#endif // GF2_MERGING_HPP
