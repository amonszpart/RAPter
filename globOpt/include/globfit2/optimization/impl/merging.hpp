#ifndef GF2_MERGING_HPP
#define GF2_MERGING_HPP

#include "pcl/console/parse.h"

#include "globfit2/optimization/merging.h"
#include "globfit2/optimization/mergingFunctors.h"
#include "globfit2/optimization/energyFunctors.h"

#include "globfit2/parameters.h"
#include "globfit2/visualization/visualization.h"
#include "globfit2/io/io.h"
#include "globfit2/processing/util.hpp"          //getPopulations()
#include "globfit2/optimization/patchDistanceFunctors.h" // RepresentativeSqrPatchPatchDistanceFunctorT
#include "globfit2/util/util.hpp"

#define CHECK(err,text) { if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << text << " returned an error! Code: " << err << std::endl; }

namespace GF2 {

namespace merging
{
    /*! \brief Dummy concept, how to use \ref processing::transformPrimitiveMap in globfit2/processing/util.hpp.
     */
    template <class _PrimitiveT, class _PointContainerT, typename _Scalar>
    struct RefitFunctor
    {
            RefitFunctor( _PointContainerT const& points, GidPidVectorMap const& populations, _Scalar scale ) : _points(points), _populations(populations), _scale(scale) {}
            int eval( _PrimitiveT& prim ) const
            {
                const int gid = prim.getTag( _PrimitiveT::GID );
                const int dir_gid = prim.getTag( _PrimitiveT::DIR_GID );
                std::cout << "[" << __func__ << "]: "
                          << "transforming primitive with GID: " << gid
                          << ", and DIR_GID: " << dir_gid
                          << std::endl;

                GidPidVectorMap::const_iterator pop_it = _populations.find( prim.getTag(_PrimitiveT::GID) );
                if ( (pop_it != _populations.end()) && (pop_it->second.size()) ) // if population exists, and has non-zero points
                {
                    std::cout << "refit from: " << prim().transpose();
                    _PrimitiveT refit;
                    processing::fitLinearPrimitive<_PrimitiveT::Dim>( /* [in,out] primitives: */ refit
                                                                    , /*              points: */ _points
                                                                    , /*               scale: */ _scale
                                                                    , /*             indices: */ &(pop_it->second)
                                                                    , /*    refit iter count: */ 2                 // fit and refit twice
                                                                    , /*    start from input: */ &prim             // use to calculate initial weights
                                                                    , /* refit position only: */ false
                                                                    , /*               debug: */ false  );
                    // save tags
                    Taggable tmp; tmp.copyTagsFrom( prim );
                    // create new primitive
                    prim = _PrimitiveT( refit.template pos(), prim.template dir() );
                    // rewrite tags
                    prim.copyTagsFrom( tmp );

                    std::cout << " to: "  << prim().transpose();
                    std::cout << " from " << _populations.at(gid).size() << " points" << std::endl;

                    return 0;
                }
                else
                    return 0;
            } //...eval()

            _PointContainerT const& _points;
            GidPidVectorMap  const& _populations;
            _Scalar                 _scale;
    };

    /*! \brief Dummy concept, how to use \ref processing::transformPrimitiveMap in globfit2/processing/util.hpp.
     */
    template <class _PrimitiveT>
    struct DirectionGroupArityFunctor
    {
            inline
            DirectionGroupArityFunctor( ) {}

            inline
            int eval(const  _PrimitiveT& prim, int )
            {
                _arities[prim.getTag( _PrimitiveT::DIR_GID )] ++;
                return 0;
            } //...eval()

            std::map<int,int> _arities;
    };
}

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
    std::vector<_Scalar>  angle_gens = { _Scalar(90.) };
    // parse params
    {
        bool valid_input = true;
        params.is3D = pcl::console::find_switch(argc,argv,"--merge3D");

        valid_input &= pcl::console::parse_argument( argc, argv, "--scale", params.scale ) >= 0;
        valid_input &= pcl::console::parse_argument( argc, argv, "--prims", prims_path   ) >= 0;

        // cloud
        pcl::console::parse_argument( argc, argv, "--cloud", cloud_path );
        valid_input &= boost::filesystem::exists( cloud_path );

        pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens );
        pcl::console::parse_argument( argc, argv, "--adopt", params.do_adopt );
        pcl::console::parse_argument( argc, argv, "--thresh-mult", params.spatial_threshold_mult );
        pcl::console::parse_argument( argc, argv, "--assoc", assoc_path );
        pcl::console::parse_argument( argc, argv, "-a", assoc_path );

        std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt --formulate\n"
                  << "\t--scale " << params.scale << "\n"
                  << "\t--prims " << prims_path << "\n"
                  << "\t--cloud " << cloud_path << "\n"
                  << "\t-a,--assoc " << assoc_path << "\n"
                  << "\t[--angle-gens "; for(int i=0;i!=angle_gens.size();++i)std::cerr<<angle_gens[i]<<",";std::cerr<<"]\n";
        std::cerr << "\t[--adopt " << params.do_adopt << "]\n"
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
    processing::appendAnglesFromGenerators( params.angles, angle_gens, true );

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
        std::cout << "starting adoptPoints" << std::endl; fflush(stdout);
        if (params.is3D)
            adoptPoints<GF2::MyPointFinitePlaneDistanceFunctor, _PointPrimitiveT, _PrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>
                    ( points, prims_map, params.scale, params.do_adopt );
        else
            adoptPoints<GF2::MyPointFiniteLineDistanceFunctor, _PointPrimitiveT, _PrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>
                    ( points, prims_map, params.scale, params.do_adopt );
    }

    // MERGE
    std::cout << "starting mergeSameDirGids" << std::endl; fflush(stdout);
    RepresentativeSqrPatchPatchDistanceFunctorT<_Scalar, SpatialPatchPatchSingleDistanceFunctorT<_Scalar> >
            patchPatchDistFunct( params.scale * params.patch_dist_limit_mult
                               , params.angle_limit
                               , params.scale
                                 , params.patch_spatial_weight );

    PrimitiveMapT prims_map_copy = prims_map;
    PrimitiveMapT out_prims,
            *in  = &prims_map_copy,
            *out = &out_prims;

    //some test here, nothing really worked.
    // The idea was to try to generate the right functor to merge either lines or planes
    //auto decideMergeFunct = params.is3D ? DecideMergePlaneFunctor() : DecideMergeLineFunctor();
    //auto decideMergeFunct = *(params.is3D ? ((void*)(&(DecideMergePlaneFunctor))) : ((void*)(&(DecideMergeLineFunctor()))));
    //
    // nico: This is ugly, but I didn't find a nice way to handle that (something like lazy type evaluation)
    if (params.is3D){
    while(
          mergeSameDirGids<_PrimitiveT, _PointPrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>
          ( *out, points, *in, params.scale, params.spatial_threshold_mult * params.scale, params.parallel_limit, patchPatchDistFunct, DecideMergePlaneFunctor() ))
    {
        PrimitiveMapT* tmp = out;
        out = in;
        in  = tmp;
    }
    }else {// 2D
        while(
              mergeSameDirGids<_PrimitiveT, _PointPrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>
              ( *out, points, *in, params.scale, params.spatial_threshold_mult * params.scale, params.parallel_limit, patchPatchDistFunct, DecideMergeLineFunctor() ))
        {
            PrimitiveMapT* tmp = out;
            out = in;
            in  = tmp;
        }
    }


    // dummy example of an iteration over all primitives
    if ( 0 ) // refit produces nan-s in the lines, not sure if it's because of the merge input, or processing::fitLinearPrimitive
    {
        typedef merging::RefitFunctor<_PrimitiveT,_PointContainerT,_Scalar> RefitFunctorT;

        // get points assigned to each patch
        GidPidVectorMap populations;
        processing::getPopulations( populations, points );
        // initialize refit functor
        RefitFunctorT refitFunctor( points, populations, params.scale );
        // refit all lines
        processing::transformPrimitivesMap< _PrimitiveT
                                          , typename PrimitiveMapT::mapped_type::iterator
                                          , RefitFunctorT
                                          > ( /* [in,out] primitives: */ *out
                                            , /* [in]        functor: */ refitFunctor );
    }

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
        io::savePrimitives   <_PrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>( *out, o_path );
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
    typedef typename _PointContainerT::value_type PointT;

    int err = EXIT_SUCCESS;

    // Loop over all points, and select orphans

    bool changed = false;

    do{

        changed = false;

        _PointPrimitiveDistanceFunctor distFunctor;

        // Populations
        GidPidVectorMap populations; // populations[gid] == std::vector<int> {pid0,pid1,...}
        if ( EXIT_SUCCESS == err )
        {
            err = processing::getPopulations( populations, points );
            CHECK( err, "getPopulations" );
        }

        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            typename _PrimitiveContainerT::const_iterator it = prims.find(points[pid].getTag( _PointPrimitiveT::GID ));

            if (    ( it == prims.end())
                    || (!containers::valueOf<_PrimitiveT>(it).size()) )
            {
                // We here have an orphean, se we need to iterate over all primitives and get the closest distance < scale
                _Scalar minDist = std::numeric_limits<_Scalar>::max();
                int minGid = -1;
                auto pos = points[pid].pos();


                // now loop over all primitives
                for ( outer_const_iterator it1 = prims.begin(); it1 != prims.end(); ++it1 )
                {
                    for ( _inner_const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 )
                    {
                        int gid = (*it2).getTag(_PrimitiveT::GID);

                        // this is very unefficient, extents must be stored !
                        // ...
                        // I'm too lazy to do it, sorry (I really apologise)
                        // ...
                        // Yes I know there is such code 10 lines below.
                        // ...
                        // Just do it and leave me alone
                        std::vector< Eigen::Matrix<_Scalar,3,1> > extrema;
                        err = it2->template getExtent<_PointPrimitiveT>
                                ( extrema
                                  , points
                                  , scale
                                  , populations[gid].size() ? &(populations[gid]) : NULL );

                        _Scalar dist = distFunctor.eval(extrema, *it2, pos);
                        if (dist < minDist){
                            minDist = dist;
                            minGid  = gid;
                        }

                    }
                }

                //std::cout << "Process orphan.... " << minDist << std::endl;


                if(minDist < scale && minDist>=0){
                    points[pid].setTag( _PointPrimitiveT::GID, minGid );
                    std::cout << "Orphan re-assigned " << pid << " " << minGid << std::endl;
                    changed = true;
                }
            }
        }
    } while (changed);

#if 0
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
#endif
} //...adoptPoints()
template <class Container,
          class PrimitiveT,
          class Population,
          class PointCloud,
          class Scalar>
inline void merge( Container&        out_primitives, // [out] Container storing merged primitives
                   PrimitiveT        l0,             // [in]  First primitive (local copy)
                   const Population& pop0,           // [in]  First primitive population (point ids)
                   PrimitiveT        l1,             // [in]  Second primitive (local copy)
                   const Population& pop1,           // [in]  Second primitive population (point ids)
                   PointCloud&       points,         // [in]  Point cloud
                   Scalar            scale,          // [in]  Working scale (for refit)
                   bool              nogeneration=false){

    int gid0 = l0.getTag(PrimitiveT::GID );
    int gid1 = l1.getTag(PrimitiveT::GID );
    int did0 = l0.getTag(PrimitiveT::DIR_GID );
    int did1 = l1.getTag(PrimitiveT::DIR_GID );

    std::cout << "merging (" << gid0 << ","
                             << did0 << ") "
//                             << l0.pos().transpose() <<  " - "
//                             << l0.dir().transpose()
                             <<  std::endl
              << " and    (" << gid1 << ","
                             << did1 << ") "
//                             << l1.pos().transpose() <<  " - "
//                             << l1.dir().transpose()
                             <<  std::endl;


    // Prepare population for refit (pop + pop1)
    Population pop = pop0;
    pop.insert(pop.end(), pop1.begin(), pop1.end());

    if (pop.size() == 0){
        cerr << "THERE IS SOMETHING WRONG HERE" << endl;
        exit(-10);
    }

    // When the direction id of the two primitives are differents,
    // the merging behavior is changing according to the arity of the
    // involved direction groups in which the two primitives belong:
    //  a. two groups arity = 1: we have no directionnal constraint on both primitives, so
    //                           we merge the two populations, refit position and direction,
    //                           store all of this in the first primitive and discard the
    //                           second one.
    //  b. one group arity = 1:  the associated primitive can be discarded, and its points merged into the
    //                           other one. A position-only fit must then be recomputed to preserve
    //                           the directionnal constraint (group arity != 1)
    //  c. no group arity = 1:   both primitives are constrained by direction. Here we duplicate both of
    //                           them, and insert them to the other direction group (position fit only).
    //
    // When the direction id of the two primitives are equals, we can merge refit only the position.
    // This lead to call case b.


    // Here are the ids that we be used to recompute the assignement after refit:
    // points with GID=originalGid will be re-assigned to newGid.
    // By default, we always transferts the points from l1 to l0 (see below)
    // but this can be changed regarding the case a-b-c.
    int originalGid = gid1;
    int newGid      = gid0;

    //std::cout << "compute Arity...." << std::endl;

    // Compute arities
    merging::DirectionGroupArityFunctor<PrimitiveT> functor;
    processing::filterPrimitives<PrimitiveT,
            typename Container::mapped_type::const_iterator > (out_primitives, functor);
    std::cout << "compute Arity....DONE" << std::endl;
    int arity0 = functor._arities[did0];
    int arity1 = functor._arities[did1];


    //std::cout << "remove previous instances" << std::endl;

    // Now we can remove the two primitives from the container
    // They we be replaced in the next part of the function, and we have the
    // local copies l0 and l1 to get their properties
    {
        typename Container::mapped_type& innerContainer0 = out_primitives.at(gid0);
        typename Container::mapped_type& innerContainer1 = out_primitives.at(gid1);
        int uid0 = l0.getTag(PrimitiveT::USER_ID1);
        int uid1 = l1.getTag(PrimitiveT::USER_ID1);

        for(typename Container::mapped_type::iterator it = innerContainer0.begin();
            it != innerContainer0.end(); it++)
            if((*it).getTag(PrimitiveT::USER_ID1) == uid0){
                innerContainer0.erase(it);
                cout << "erase me" << endl;
                break;
            }
        if (innerContainer0.size() == 0)  out_primitives.erase(gid0);
        for(typename Container::mapped_type::iterator it = innerContainer1.begin();
            it != innerContainer1.end(); it++)
            if((*it).getTag(PrimitiveT::USER_ID1) == uid1){
                innerContainer1.erase(it);
                cout << "erase me" << endl;
                break;
            }
        if (innerContainer1.size() == 0) out_primitives.erase(gid1);
    }


    // temporary array containing generated primitives
    std::vector<PrimitiveT> primToAdd;

    if (arity0 == 1 && arity1 == 1){
        // two free patches, that have no constraints on their direction
        // here the plan is
        std::cout << "Case A" << std::endl;

        PrimitiveT mergedPrim;
        processing::fitLinearPrimitive<PrimitiveT::Dim>( /* [in,out] primitives: */ mergedPrim
                                                       , /*              points: */ points
                                                       , /*               scale: */ scale
                                                       , /*             indices: */ &pop
                                                       , /*    refit iter count: */ 5);               // fit and refit twice

        // by default we copy from l0
        mergedPrim.copyTagsFrom(l0);
        primToAdd.push_back(mergedPrim);

    }else /* if (arity0 != 1 && arity1 != 1)*/{
        // two constrained patches, we need to generate new primitives.
        std::cout << "Case C: " << arity0 << " - " << arity1 << std::endl;

        // Compute population centroid.
        auto centroid = processing::getCentroid<Scalar>(points, &pop);

        // Copy from l0
        PrimitiveT mergedPrim (centroid, l0.dir());
        mergedPrim.copyTagsFrom(l0);
        primToAdd.push_back(mergedPrim);

        // Copy from l1
        mergedPrim = PrimitiveT(centroid, l1.dir());
        mergedPrim.copyTagsFrom(l1);
        // set same group id as l0 (same population)
        mergedPrim.setTag(PrimitiveT::GID, l0.getTag(PrimitiveT::GID) );
        primToAdd.push_back(mergedPrim);

    }/*else {
        // one patch is free, the other constrained.
        std::cout << "Case B" << std::endl;

        PrimitiveT *sourcePrim = arity0 != 1 ? & l0 : & l1;

        // refit position
        PrimitiveT mergedPrim (processing::getCentroid<Scalar>(points, &pop), sourcePrim->dir());
        mergedPrim.copyTagsFrom(*sourcePrim);
        primToAdd.push_back(mergedPrim);

        // this is a constrained fit so we may need to invert the arity transfer direction
        if (arity0 == 1){
            originalGid = gid0;
            newGid      = gid1;
        }
    }*/

    // add generated primitives
    for(typename std::vector<PrimitiveT>::const_iterator it = primToAdd.begin(); it != primToAdd.end(); it++){
        containers::add<PrimitiveT>( out_primitives, (*it).getTag(PrimitiveT::GID ), (*it) );
        std::cout<< "Add (" << (*it).getTag(PrimitiveT::GID )
                 << ","     << (*it).getTag(PrimitiveT::DIR_GID ) << ") : "
//                 << (*it).pos().transpose() <<  " - "
//                 << (*it).dir().transpose() <<  " - "
                 << std::endl;
    }

    //std::cout << "recompute assignment: Point " << originalGid << " will now be " << newGid << std::endl;

    // Recompute assignement, from originalGid to newGid
    typedef typename PointCloud::value_type PointT;
    for ( size_t pid = 0; pid != points.size(); ++pid )
    {
        if ( points[pid].getTag( PointT::GID ) == originalGid){
            points[pid].setTag(  PointT::GID, newGid );
        }
    }

    // for now, do nothing except: the two primitives are kept

    //containers::add<PrimitiveT>( out_primitives, l0.getTag(PrimitiveT::GID ), l0 );
    //containers::add<PrimitiveT>( out_primitives, l1.getTag(PrimitiveT::GID ), l1 );
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
         , class    _PatchPatchDistanceFunctorT
         , class    _PrimitiveDecideMergeFunctorT >
int Merging::mergeSameDirGids( _PrimitiveContainerT             & out_primitives
                             , _PointContainerT                 & points
                             , _PrimitiveContainerT        /*const&*/ primitives
                             , _Scalar                     const  scale
                             , _Scalar                     const  spatial_threshold
                             , _Scalar                     const  parallel_limit
                             , _PatchPatchDistanceFunctorT const& patchPatchDistFunct
                             , _PrimitiveDecideMergeFunctorT const& primitiveDecideMergeFunct  )
{
    typedef typename _PrimitiveContainerT::const_iterator      outer_const_iterator;
    typedef           std::vector<Eigen::Matrix<_Scalar,3,1> > ExtremaT;
    typedef           std::map   < int, ExtremaT>              LidExtremaT;
    typedef           std::map   < int, LidExtremaT >          GidLidExtremaT;
    typedef           std::pair  < int, // map key
                                   int> // linear index in the array associated to the key
                      GidLid;



    // Store the primitives that have been matched and must be ignored
    // First  (key)   = candidate,
    // Second (value) = reference.
    typedef std::set< GidLid > IgnoreListT;
    IgnoreListT ignoreList;


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
                                  (outer_it != primitives.end()); // we now handle error
                                 ++outer_it )
        {
            int gid  = -2; // (-1 is "unset")
            int lid = 0; // linear index of primitive in container (to keep track)
            // for all directions
            for ( _inner_const_iterator inner_it  = containers::valueOf<_PrimitiveT>(outer_it).begin();
                                       (inner_it != containers::valueOf<_PrimitiveT>(outer_it).end());// we now handle error
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
                                                           , points
                                                           , scale
                                                           , populations[gid].size() ? &(populations[gid]) : NULL );

                if(err != EXIT_SUCCESS){
                    std::cerr << "Issue when computing extent of ("
                              << primitives.at(gid).at(lid).getTag(_PrimitiveT::GID )     << ","
                              << primitives.at(gid).at(lid).getTag(_PrimitiveT::DIR_GID ) << ")"
                              << std::endl
                              << "Ignored later... " << std::endl;

                    ignoreList.insert(GidLid (gid, lid));
                }
            } //...for primitives
        } //...for patches

        CHECK( err, "calcExtrema" )
    } //...getExtrema




    typedef typename GidLidExtremaT::const_iterator GidIt;
    typedef typename LidExtremaT::const_iterator    PrimIt;

    // add a local uid attribute to primitives to be able to recognize them later
    int uid = 0;
    for ( GidIt gid_it = extrema.cbegin(); gid_it != extrema.cend(); ++gid_it ){
        const int gid0 = gid_it->first;
        for ( PrimIt prim_it = gid_it->second.cbegin(); prim_it != gid_it->second.cend(); ++prim_it, uid++ ){
            const int lid0 = std::distance<typename LidExtremaT::const_iterator>( gid_it ->second.cbegin(), prim_it );
            primitives.at(gid0).at(lid0).setTag(_PrimitiveT::USER_ID1, uid);
        }
    }

//    // debug
//    cout << "[in]:  " << primitives.size() << endl;
//    for ( outer_const_iterator outer_it  = primitives.begin();
//          (outer_it != primitives.end());
//          ++outer_it )
//        for(auto innerIt = (*outer_it).second.begin(); innerIt != (*outer_it).second.end(); innerIt ++)
//            cout << (*innerIt).getTag(_PrimitiveT::GID ) << " - "
//                 << (*innerIt).getTag(_PrimitiveT::DIR_GID ) << endl;

    // Here are two first loops to iterate over the map entries, and for each of them over the
    // linear array containing primitives (we call them ref. primitive in the following).
    //
    // For each of these ref. primitives, we check for merging candidates. They can be in:
    //   - the same map entry, after the ref. primitive in the linear Primitive array
    //   - the next map entries, any index
    //
    // For each couple ref/candidate, we check if we can merge. If yes, we do it and then invalidate
    // both the ref and the candidate to prevent to merge them with other primitives. Indeed, the
    // merging process can potentially remove the primitives, or at least change their properties.
    // Primitives are invalidated by storing them as key in the ignoreList structure.
    //
    // After processing, the ignoreList contains any candidate described by a GidLid gl:
    //   ignoreList [gl];
    //
    // The output buffer is initialized with the input. All merging operations will remove
    // old primitives and replace them by merged one.
    out_primitives = primitives;

    // Reference traversal
    for ( GidIt gid_it = extrema.cbegin(); gid_it != extrema.cend(); ++gid_it )
    {
        const int gid0 = gid_it->first;

        for ( PrimIt prim_it = gid_it->second.cbegin(); prim_it != gid_it->second.cend(); ++prim_it )
        {
            // linear id of the reference
            const int lid0 = std::distance<typename LidExtremaT::const_iterator>( gid_it ->second.cbegin(), prim_it );

            // Reference primitive key, used to check if it has not been merged previously
            GidLid refGidLid (gid0, lid0);

            // check if this primitives has not been merged previously
            if (ignoreList.find(refGidLid) != ignoreList.end()) continue;

            // reference primitive
            const _PrimitiveT& prim0 = primitives.at(gid0).at(lid0);

            // is the reference valid, i.e. no merge occured
            bool is0Valid = true;


            // Candidates traversal
            // This loop and the nested one check if the ref primitive is still valid for comparison,
            // i.e. no merge occured.
            // This construction prevent to use two successive break commands (this loop and the nested one)
            for ( GidIt gid_it1 = gid_it; is0Valid && gid_it1 != extrema.cend(); ++gid_it1 )
            {

                int gid1 = gid_it1->first;

                // Given the ref. primitive, we check for merging candidates. They can be in:
                //   - the same map entry, after the ref. primitive in the linear Primitive array
                //   - the next map entries, any index
                for ( PrimIt prim_it1 = (gid_it1 == gid_it) ? ++PrimIt(prim_it) :       // for the same map entry, after the ref. primitive
                                                              gid_it1->second.cbegin(); // for the next map entries, any index
                      is0Valid && prim_it1 != gid_it1->second.cend();
                      ++prim_it1 )
                {

                    const int lid1 = std::distance<typename LidExtremaT::const_iterator>( gid_it1->second.cbegin(), prim_it1 );

                    // resume of the candidate primitive, used to record the matching operation for later and
                    // detect previous merge operations
                    GidLid candGidLid (gid1, lid1);

                    // Here we don't need to define
                    // bool is1Valid,
                    // calling continue is sufficient to jump to the next primitive after and merge,
                    // plus here check that a previous merge has not been recorded
                    if (ignoreList.find(candGidLid) != ignoreList.end()) continue;

                    const _PrimitiveT& prim1 = primitives.at(gid1).at(lid1);

                    if (primitiveDecideMergeFunct.eval( prim_it->second,  // extrema 0
                                                  prim0,            // prim 0
                                                  prim_it1->second, // extrema 1
                                                  prim1,            // prim 1
                                                  scale))
                    {
                        //std::cout << " YES" << std::endl;

                        // record this to detect unmerged primitives later and invalidate both primitives
                        ignoreList.insert(candGidLid);
                        ignoreList.insert(refGidLid);

                        merge( out_primitives,     // [out] Container storing merged primitives
                               prim0,              // [in]  First primitive (can be invalidated during the call)
                               populations[gid0],  // [in]  First primitive population (point ids)
                               prim1,              // [in]  Second primitive
                               populations[gid1],  // [in]  Second primitive population (point ids)
                               points,             // [in]  Point cloud
                               scale);             // [in]  Working scale (for refit)

                        is0Valid = false;

                        continue;  // jump to the next couple ref/candidate primitive
                    }
                    //else
                    //    std::cout << " NO" << std::endl;
                }
            }
        }
    }

    typedef typename _PrimitiveContainerT::mapped_type::iterator inner_iterator;

    // we can now remove primitives that are not assigned to any points
    processing::eraseNonAssignedPrimitives<_PrimitiveT, inner_iterator>(out_primitives, points);

//    cout << "[out]: " << out_primitives.size() << endl;
//    for ( outer_const_iterator outer_it  = out_primitives.begin();
//          (outer_it != out_primitives.end()) ;
//          ++outer_it )
//        for(auto innerIt = (*outer_it).second.begin(); innerIt != (*outer_it).second.end(); innerIt ++)
//            cout << (*innerIt).getTag(_PrimitiveT::GID ) << " - "
//                 << (*innerIt).getTag(_PrimitiveT::DIR_GID ) << endl;


    // here we detect if merge have occured and changed the shape of the primitive container
    // Indeed, sometime the merge can generate 2 primitives that can be merged the next times,
    // leading to an infinite loop.
    // So we check that the primitive arrays of each gid are different
    if(primitives.size() != out_primitives.size())
        return true; // we need another iteration

    typename _PrimitiveContainerT::iterator outer_itIn  =     primitives.begin();
    typename _PrimitiveContainerT::iterator outer_itOut = out_primitives.begin();


    auto cmp_primitive = [](_PrimitiveT const& l0, _PrimitiveT const& l1)
    {
        if (l0.getTag(_PrimitiveT::GID ) == l1.getTag(_PrimitiveT::GID ))
            return l0.getTag(_PrimitiveT::DIR_GID ) < l1.getTag(_PrimitiveT::DIR_GID );
        return l0.getTag(_PrimitiveT::GID ) < l1.getTag(_PrimitiveT::GID );
    };

    for ( ;
          (outer_itIn != primitives.end()) ; // don't need to check for out, we know they have the same size
          ++outer_itIn, ++outer_itOut ){
        if ( (*outer_itIn ).second.size() != (*outer_itOut ).second.size() )
            return true; // primitive inserted or deleted

        std::sort ((*outer_itIn ).second.begin(), (*outer_itIn ).second.end(), cmp_primitive);
        std::sort ((*outer_itOut).second.begin(), (*outer_itOut).second.end(), cmp_primitive);

        inner_iterator inner_itIn  = (*outer_itIn ).second.begin();
        inner_iterator inner_itOut = (*outer_itOut).second.begin();

        for( ; inner_itIn != (*outer_itIn ).second.end(); ++inner_itIn, ++inner_itOut){
            _PrimitiveT& l0 = *inner_itIn;
            _PrimitiveT& l1 = *inner_itOut;

            if (l0.getTag(_PrimitiveT::DIR_GID ) != l1.getTag(_PrimitiveT::DIR_GID ) ||
                l0.getTag(_PrimitiveT::GID ) != l1.getTag(_PrimitiveT::GID ) )
                return true;
        }
    }

    // if we reach that point, that means that input/output are similar
    return false;

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
} //...Merging::mergeSameDirGids()

} //...namespace GF2

#undef CHECK

#endif // GF2_MERGING_HPP
