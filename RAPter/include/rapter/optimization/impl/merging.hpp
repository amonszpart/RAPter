#ifndef RAPTER_MERGING_HPP
#define RAPTER_MERGING_HPP

#include <iostream>

#include "rapter/util/parse.h"

#include "rapter/optimization/merging.h"
#include "rapter/optimization/mergingFunctors.h"
#include "rapter/optimization/energyFunctors.h"

#include "rapter/parameters.h"
//#include "rapter/visualization/visualization.h"
#include "rapter/io/io.h"
#include "rapter/processing/util.hpp"          //getPopulations()
#include "rapter/processing/impl/angleUtil.hpp" // appendAngles...
#include "rapter/optimization/patchDistanceFunctors.h" // RepresentativeSqrPatchPatchDistanceFunctorT
#include "rapter/util/util.hpp"
#include "omp.h"

#define CHECK(err,text) { if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << text << " returned an error! Code: " << err << std::endl; }

namespace rapter {

namespace merging
{
    /*! \brief Dummy concept, how to use \ref processing::transformPrimitiveMap in rapter/processing/util.hpp.
     */
    template <class _PrimitiveT, class _PointContainerT, typename _Scalar>
    struct RefitFunctor
    {
            RefitFunctor( _PointContainerT const& points, GidPidVectorMap const& populations, _Scalar scale ) : _points(points), _populations(populations), _scale(scale) {}
            int eval( _PrimitiveT& prim ) const
            {
                const GidT gid = prim.getTag( _PrimitiveT::TAGS::GID );
                const DidT dir_gid = prim.getTag( _PrimitiveT::TAGS::DIR_GID );
                std::cout << "[" << __func__ << "]: "
                          << "transforming primitive with GID: " << gid
                          << ", and DIR_GID: " << dir_gid
                          << std::endl;

                GidPidVectorMap::const_iterator pop_it = _populations.find( prim.getTag(_PrimitiveT::TAGS::GID) );
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
                    typename _PrimitiveT::TaggableT tmp; tmp.copyTagsFrom( prim );
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

    /*! \brief Dummy concept, how to use \ref processing::transformPrimitiveMap in rapter/processing/util.hpp.
     */
    template <class _PrimitiveT>
    struct DirectionGroupArityFunctor
    {
            inline
            DirectionGroupArityFunctor( ) {}

            inline
            LidT eval(const  _PrimitiveT& prim, LidT )
            {
                _arities[prim.getTag( _PrimitiveT::TAGS::DIR_GID )] ++;
                return 0;
            } //...eval()

            std::map<DidT,LidT> _arities;
    };

    /*! \brief Merge two primitives, after the decision has been taken.
     *  \param[in] max_dir_gid Holds the current maximum dir_gid that's taken. A new direction should get \p max_dir_gid +1.
     */
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
                       DidT            & max_dir_gid
                     , PidT            & maxUid4         // [in/out] Keeps track of new unique ids, that  live outside the big while loop
                     , bool              nogeneration = false
                     )
    {

        GidT gid0 = l0.getTag(PrimitiveT::TAGS::GID );
        GidT gid1 = l1.getTag(PrimitiveT::TAGS::GID );
        DidT did0 = l0.getTag(PrimitiveT::TAGS::DIR_GID );
        DidT did1 = l1.getTag(PrimitiveT::TAGS::DIR_GID );

        // Prepare population for refit (pop + pop1)
        Population pop = pop0;
        pop.insert(pop.end(), pop1.begin(), pop1.end());

        if (pop.size() == 0){
            std::cerr << "[" << __func__ << "]: " << "THERE IS SOMETHING WRONG HERE: pop.size() == 0" << std::endl;
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
        GidT originalGid = gid1;
        GidT newGid      = gid0;

        //std::cout << "compute Arity...." << std::endl;

        // Compute arities
        merging::DirectionGroupArityFunctor<PrimitiveT> functor;
        processing::filterPrimitives<PrimitiveT,
                typename Container::mapped_type/*::const_iterator*/ > (out_primitives, functor);
    //    std::cout << "compute Arity....DONE" << std::endl;
        int arity0 = functor._arities[did0];
        int arity1 = functor._arities[did1];


        //std::cout << "remove previous instances" << std::endl;

        // Now we can remove the two primitives from the container
        // They we be replaced in the next part of the function, and we have the
        // local copies l0 and l1 to get their properties
        {
            typename Container::mapped_type& innerContainer0 = out_primitives.at(gid0);
            typename Container::mapped_type& innerContainer1 = out_primitives.at(gid1);
            UidT uid0 = l0.getTag(PrimitiveT::USER_ID1);
            UidT uid1 = l1.getTag(PrimitiveT::USER_ID1);

            for(typename Container::mapped_type::iterator it = innerContainer0.begin();
                it != innerContainer0.end(); it++)
                if((*it).getTag(PrimitiveT::USER_ID1) == uid0){
                    innerContainer0.erase(it);
    //                cout << "erase me" << endl;
                    break;
                }
            if (innerContainer0.size() == 0)  out_primitives.erase(gid0);
            for(typename Container::mapped_type::iterator it = innerContainer1.begin();
                it != innerContainer1.end(); it++)
                if((*it).getTag(PrimitiveT::USER_ID1) == uid1){
                    innerContainer1.erase(it);
    //                cout << "erase me" << endl;
                    break;
                }
            if (innerContainer1.size() == 0) out_primitives.erase(gid1);
        }


        // temporary array containing generated primitives
        std::vector<PrimitiveT> primToAdd;

        if (arity0 == 1 && arity1 == 1){
            // two free patches, that have no constraints on their direction
            // here the plan is
    //        std::cout << "Case A" << std::endl;

            PrimitiveT mergedPrim;
            processing::fitLinearPrimitive<PrimitiveT::Dim>( /* [in,out] primitives: */ mergedPrim
                                                           , /*              points: */ points
                                                           , /*               scale: */ scale
                                                           , /*             indices: */ &pop
                                                           , /*    refit iter count: */ 5);               // fit and refit twice

            // by default we copy from l0
            mergedPrim.copyTagsFrom(l0);
            // change direction ID to the next free one, since we changed the direction with refit
            mergedPrim.setTag( PrimitiveT::TAGS::DIR_GID, ++max_dir_gid );
            // save
            primToAdd.push_back(mergedPrim);

        }else /* if (arity0 != 1 && arity1 != 1)*/{
            // two constrained patches, we need to generate new primitives.
    //        std::cout << "Case C: " << arity0 << " - " << arity1 << std::endl;

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
            mergedPrim.setTag(PrimitiveT::TAGS::GID, l0.getTag(PrimitiveT::TAGS::GID) );
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
        for(typename std::vector<PrimitiveT>::iterator it = primToAdd.begin(); it != primToAdd.end(); it++)
        {
            (*it).template setExtentOutdated();
            (*it).setTag( PrimitiveT::USER_TAGS::USER_ID4, ++maxUid4 );
            containers::add<PrimitiveT>( out_primitives, (*it).getTag(PrimitiveT::TAGS::GID ), (*it) );
    //        std::cout<< "Add (" << (*it).getTag(PrimitiveT::TAGS::GID )
    //                 << ","     << (*it).getTag(PrimitiveT::TAGS::DIR_GID ) << ") : "
    ////                 << (*it).pos().transpose() <<  " - "
    ////                 << (*it).dir().transpose() <<  " - "
    //                 << std::endl;
        }

        //std::cout << "recompute assignment: Point " << originalGid << " will now be " << newGid << std::endl;

        // Recompute assignement, from originalGid to newGid
        typedef typename PointCloud::value_type PointT;
        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            if ( points[pid].getTag( PointT::TAGS::GID ) == originalGid){
                points[pid].setTag(  PointT::TAGS::GID, newGid );
            }
        }

        // for now, do nothing except: the two primitives are kept

        //containers::add<PrimitiveT>( out_primitives, l0.getTag(PrimitiveT::TAGS::GID ), l0 );
        //containers::add<PrimitiveT>( out_primitives, l1.getTag(PrimitiveT::TAGS::GID ), l1 );
    }


    template <typename _ULongT = ULidT>
    class ComparedSet : public std::set< std::pair<_ULongT, _ULongT> >
    {
        public:
            typedef std::pair<_ULongT, _ULongT> ElementT;
            typedef std::set < ElementT       > ParentT;

            /*! \brief Overload constructor to initialize max id field.
             */
            ComparedSet() : ParentT(), _maxId( _ULongT(0) ), _hits(0) {}

            /*! \brief Overload set insert to keep track of maximum node id in the edges list.
             */
            std::pair<typename ParentT::iterator, bool>
            insert( typename ParentT::value_type& __x )
            {
                // this is just for safety, _maxId is inited, changed from the outside
                if ( __x.first  > _maxId ) _maxId = __x.first;
                if ( __x.second > _maxId ) _maxId = __x.second;

                return ParentT::insert( __x );
            }

            /*! \brief Erases any member pair, that has x.first or x.second in it.
             */
            void eraseAny( typename ParentT::value_type& x )
            {
                for ( auto it = this->begin(); it != this->end(); )
                    if (    (it->first == x.first) || (it->second == x.first)
                         || (it->first== x.second) || (it->second == x.second)
                       )
                        this->erase( it++ );
                    else
                        ++it;
            }

            void
            clear() _GLIBCXX_NOEXCEPT
            {
                ParentT::clear();
                _maxId = 0;
            }

            /*! \brief Get maximum vertex id, this is changed from outside, not the "seen type"
             */
            inline _ULongT& getMaxId() { return _maxId; }
            inline void incHits() { _hits++; }
            inline ULidT getHits() { return _hits; }

        protected:
            _ULongT _maxId;
            ULidT _hits;
    }; //...ComparedSet


    /*! \brief  Function that performs the merge as a single operation, called once from MergCli.
     *          When splitting the scene, we can call this for each split.
     * \note    Added by Aron on 12/1/2015, 16:30
     */
    template <
           class    _PointPrimitiveT
         , class    _PrimitiveT
         , class    _PointContainerT
         , class    _PrimitiveMapT
         , class    _ParamsT
         >
    void inline iterativeMerge( _PrimitiveMapT        & out_prims
                              , _PointContainerT      & points
                              , _PrimitiveMapT   const& prims_map
                              , _ParamsT         const& params
                              )
    {
        typedef typename _PrimitiveT::Scalar Scalar;
        //typedef std::vector<_PrimitiveT>    PatchT;
        //typedef std::map   < GidT, PatchT >  PrimitiveMapT;

        //_PrimitiveContainerT prims;

        // MERGE
        std::cout << "starting mergeSameDirGids" << std::endl; fflush(stdout);
        RepresentativeSqrPatchPatchDistanceFunctorT<Scalar, SpatialPatchPatchSingleDistanceFunctorT<Scalar> >
                patchPatchDistFunct( params.scale * params.patch_dist_limit_mult
                                   , params.angle_limit
                                   , params.scale
                                   , params.patch_spatial_weight
                                   );

        _PrimitiveMapT       prims_map_copy = prims_map;
        ComparedSet<PidT>    compared; // contains USER_ID4 tags, that were checked, and not merged
        for ( typename _PrimitiveMapT::Iterator it(prims_map_copy); it.hasNext(); it.step() )
        {
            it->setTag( _PrimitiveT::USER_TAGS::USER_ID4, compared.getMaxId()++ );
        }

        _PrimitiveMapT *in  = &prims_map_copy,
                       *out = &out_prims;

        //some test here, nothing really worked.
        // The idea was to try to generate the right functor to merge either lines or planes
        //auto decideMergeFunct = params.is3D ? DecideMergePlaneFunctor() : DecideMergeLineFunctor();
        //auto decideMergeFunct = *(params.is3D ? ((void*)(&(DecideMergePlaneFunctor))) : ((void*)(&(DecideMergeLineFunctor()))));
        //
        // nico: This is ugly, but I didn't find a nice way to handle that (something like lazy type evaluation)
        bool preserveSmallPatches = !params.do_adopt;

        if ( params.is3D )
        {
            while (
                    Merging::mergeSameDirGids<_PrimitiveT, _PointPrimitiveT, typename _PrimitiveMapT::mapped_type/*::const_iterator*/>
                    ( *out, points, *in, params.scale, params.spatial_threshold_mult * params.scale, params.parallel_limit
                    , patchPatchDistFunct, DecideMergePlaneFunctor(), preserveSmallPatches, compared )
                  )
            {
                _PrimitiveMapT* tmp = out;
                out = in;
                in  = tmp;
            } //...while
        } //...3D
        else
        {// 2D
            while(
                  Merging::mergeSameDirGids<_PrimitiveT, _PointPrimitiveT, typename _PrimitiveMapT::mapped_type/*::const_iterator*/>
                  ( *out, points, *in, params.scale, params.spatial_threshold_mult * params.scale, params.parallel_limit, patchPatchDistFunct, DecideMergeLineFunctor(), preserveSmallPatches, compared ))
            {
                _PrimitiveMapT* tmp = out;
                out = in;
                in  = tmp;
            } //...while
        } //...else2D

        // if last iteration output to prims_map_copy, make sure to copy to out_prims.
        if ( in != &out_prims )
        {
            std::cout << "[" << __func__ << "]: " <<  "PFFF...swapping in/out..." << std::endl;
            out_prims = *in;
        }

        std::cout << "[" << __func__ << "]: " << "UID4 HITS: " << compared.getHits() << std::endl;

    } //...doMerge cand

    template <class _PrimitiveMapT, class _PointContainerT>
    struct MergePartition : std::pair<_PrimitiveMapT, _PointContainerT>
    {
        inline _PrimitiveMapT const& getPrimitives() const { return this->first; }
        inline _PointContainerT const& getPoints() const { return this->second; }
        inline _PrimitiveMapT      & getPrimitives() { return this->first; }
        inline _PointContainerT      & getPoints()  { return this->second; }
    };

    template <class _PrimitiveMapT, class _PointContainerT>
    struct MergePartitions : std::vector< MergePartition<_PrimitiveMapT, _PointContainerT> >
    {
        typedef MergePartition<_PrimitiveMapT, _PointContainerT> ElementT;
    };

    template <class _PrimitiveMapT, class _PointContainerT, class _PartitionT, class _MergeParamsT>
    void inline partition( _PartitionT           & outPartition
                         , _PrimitiveMapT   const& prims
                         , _PointContainerT const& points
                         , _MergeParamsT    const& params
                         , size_t           const  sizeLimit
                         , int              const  splitCount = 8
                         , int              const  level = 0
                         )
    {
        typedef typename _PointContainerT::value_type            PointPrimitiveT;
        typedef typename _PrimitiveMapT::mapped_type::value_type PrimitiveT;

        // (1) recursion exit condition
        if ( prims.size() < sizeLimit )
        {
            // CALL merge!
            outPartition.getPoints() = points;
            merging::iterativeMerge<PointPrimitiveT,PrimitiveT, _PointContainerT>
                    ( /* out: */ outPartition.getPrimitives(), outPartition.getPoints()
                    , /*  in: */ prims, params );
            return;
        }

        // (2.1) split primitives
        std::vector<_PartitionT> splits; splits.resize( splitCount );
        size_t       partSize     = float(prims.size()) / float(splitCount); // targeted primcount in each partition
        size_t       splitPrimCnt = 0,                                       // current primcount in current split
                     splitId      = 0;                                       // id of current split
        std::map<GidT, LidT> gidsParts;                                      // <Gid, splitId> (helps sort the points)
        // for each input primitive
        for ( typename _PrimitiveMapT::const_iterator it = prims.begin(); it != prims.end(); ++it )
        {
            // first: GID
            // second: vector<Primitive>

            // add it to the current split
            for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2 )
            {
                containers::add( splits[ splitId ].getPrimitives(), /* gid: */ it->first, /* prim: */ *it2 );
                ++splitPrimCnt;
            }
            if ( gidsParts.find( it->first ) != gidsParts.end() )
                std::cerr << "[" << __func__ << "]: " << "gid " << it->first << " already added" << std::endl;
            // record membership for points
            gidsParts[ it->first ] = splitId;
            if ( (splitPrimCnt >= partSize) && (splitId+1 < splits.size()) )
            {
                ++splitId;
                splitPrimCnt = 0;
            }
        } //...for all input primitives

        // (2.2) split points
        for ( size_t i = 0; i < points.size(); ++i )
        {
            // get GID
            const GidT pointGid = points[i].getTag( PointPrimitiveT::TAGS::GID );
            // get split for GID
            auto it = gidsParts.find( pointGid );
            if ( it != gidsParts.end() )    splits[ it->second ].getPoints().push_back( points[i] );
            else                            splits[ 0          ].getPoints().push_back( points[i] ); // if point unassigned
        } //...split points

        // debug
        {
            std::cout << "splits: ";
            for ( size_t i = 0; i != splits.size(); ++i )
                std::cout << "<" << splits[i].getPrimitives().size() << "," << splits[i].getPoints().size() << ">, ";
            std::cout << std::endl;
        }

        // (3) recurse
        std::vector<_PartitionT> processedParts; processedParts.resize( splits.size() );
        #pragma omp parallel for if(splitCount > 2)
        for ( size_t i = 0; i < splits.size(); ++i )
        {
            partition( /* out: */ processedParts[i]
                     , /*  in: */ splits[i].getPrimitives(), splits[i].getPoints()
                     , params, sizeLimit, 2, level + 1 );
        }

        // (4.1) gather to <gatheredPrimitives,outPartition.getPoints>
        _PrimitiveMapT gatheredPrimitives;
        for ( size_t i = 0; i < processedParts.size(); ++i )
        {
            gatheredPrimitives.insert( processedParts[i].getPrimitives().begin()
                                               , processedParts[i].getPrimitives().end()   );
            outPartition.getPoints().insert( outPartition.getPoints().end()
                                             , processedParts[i].getPoints().begin()
                                             , processedParts[i].getPoints().end()  );
        } //...gather

//#if 0
        // (4.2) merge <gatheredPrimitives(),outPartition.getPoints()> to <outPartition.getPrimitives(),outPartition.getPoints()>
        if ( !level )
        {
            char o_path[2048]; sprintf( o_path, "unMergedLvl%02d", level );
            io::savePrimitives<PrimitiveT, typename _PrimitiveMapT::mapped_type::const_iterator>( gatheredPrimitives
                                                                                               , std::string(o_path) + ".csv" );
            std::cout << "wrote " << o_path << std::endl;

            char oPointsPath[2048];
            sprintf( oPointsPath, "points_%s.csv", o_path );
            io::writeAssociations<PointPrimitiveT>( outPartition.getPoints(), oPointsPath );
            std::cout << "wrote " << oPointsPath << std::endl;
        }

        {
//#else
            std::cout << "[" << __func__ << "]: " << "gathering level " << level << std::endl;
            merging::iterativeMerge<PointPrimitiveT,PrimitiveT, _PointContainerT>
                    ( /* out: */ outPartition.getPrimitives(), outPartition.getPoints()
                    , /*  in: */ gatheredPrimitives, params );
            std::cout << "[" << __func__ << "]: " << "gathering level " << level << " finished" << std::endl;
        }

        std::cout << "[" << __func__ << "]: " << "primcount: " << prims.size() << " -> " << outPartition.getPrimitives().size() << std::endl;
    } //...partition
} //...merging

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
    AnglesT  angle_gens( {AnglesT::Scalar(90.)} );
    size_t sizeLimit = 0; // if >0, a non-spatial recursive partitioning will happen

    // parse params
    {
        bool valid_input = true;
        params.is3D = pcl::console::find_switch(argc,argv,"--merge3D");

        valid_input &= pcl::console::parse_argument( argc, argv, "--scale", params.scale ) >= 0;
        valid_input &= pcl::console::parse_argument( argc, argv, "--prims", prims_path   ) >= 0;

        // cloud
        pcl::console::parse_argument( argc, argv, "--cloud", cloud_path );
        valid_input &= boost::filesystem::exists( cloud_path );

        if ( pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens ) < 0 )
        {
            std::cerr << "we need comma-separated anglegens explicitly!" << std::endl;
            throw std::runtime_error("no anglegens provided");
        }
        rapter::console::parse_argument( argc, argv, "--adopt", params.do_adopt );
        rapter::console::parse_argument( argc, argv, "--thresh-mult", params.spatial_threshold_mult );
        rapter::console::parse_argument( argc, argv, "--assoc", assoc_path );
        rapter::console::parse_argument( argc, argv, "-a", assoc_path );
        rapter::console::parse_argument( argc, argv, "--patch-pop-limit", params.patch_population_limit );

        rapter::console::parse_argument( argc, argv, "--partition", sizeLimit );

        if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale, --prims are compulsory, --cloud needs to exist" << std::endl;

            std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt --formulate\n"
                      << "\t--scale " << params.scale << "\n"
                      << "\t--prims " << prims_path << "\n"
                      << "\t--cloud " << cloud_path << "\n"
                      << "\t-a,--assoc " << assoc_path << "\n"
                      << "\t[--angle-gens "; for(size_t i=0;i!=angle_gens.size();++i)std::cerr<<angle_gens[i]<<",";std::cerr<<"]\n";
            std::cerr << "\t[--adopt " << params.do_adopt << "]\n"
                      << "\t[--patch-pop-limit " << params.patch_population_limit << "]\n"
                      << "\t[--thresh-mult " << params.spatial_threshold_mult << "]\n"
                      << "\t[--no-paral]\n"
                      << "\t[--partition " << sizeLimit << "\t split scene into chunks ]\n"
                      << std::endl;

            return EXIT_FAILURE;
        }
    } // ... parse params

    // read primitives
    //typedef std::vector<_PrimitiveT>    PatchT;
    //typedef std::map   < GidT, PatchT >  PrimitiveMapT;
    typedef containers::PrimitiveContainer<_PrimitiveT> PrimitiveMapT;
    typedef typename PrimitiveMapT::InnerContainerT     PatchT;
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
    bool no_paral = pcl::console::find_switch(argc,argv,"--no_paral");
    angles::appendAnglesFromGenerators( params.angles, angle_gens, no_paral, true );

    // associations
    std::vector<std::pair<PidT,LidT> > points_primitives;
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
        points[i].setTag( _PointPrimitiveT::TAGS::GID, points_primitives[i].first );

        // error check 2
        //if ( points[i].getTag(_PointPrimitiveT::TAGS::GID) == -1 )
        //    std::cerr << "[" << __func__ << "]: " << "point assigned to patch with id -1" << std::endl;
    }

    //____________________________WORK____________________________

    // ADOPT
    if ( params.do_adopt )
    {
        std::cout << "starting adoptPoints" << std::endl; fflush(stdout);
        if (params.is3D)
            adoptPoints<rapter::MyPointFinitePlaneDistanceFunctor, _PointPrimitiveT, _PrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>
                    ( points, prims_map, params.scale, params.do_adopt, params.patch_population_limit );
        else
            adoptPoints<rapter::MyPointFiniteLineDistanceFunctor, _PointPrimitiveT, _PrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>
                    ( points, prims_map, params.scale, params.do_adopt, params.patch_population_limit );
    }



    PrimitiveMapT out_prims;
    if ( (sizeLimit == 0) || (prims_map.size() <= sizeLimit) ) // Assumes one primitive in each gid...
        merging::iterativeMerge<_PointPrimitiveT,_PrimitiveT, _PointContainerT>
                ( /* out: */ out_prims, points, /* in: */ prims_map, params );
    else
    {
        merging::MergePartition<PrimitiveMapT, _PointContainerT> outPartition;
        partition( outPartition, prims_map, points, params, sizeLimit );
        out_prims = outPartition.getPrimitives();
        points = outPartition.getPoints();
    }

    // SAVE
    std::string o_path;
    int         iteration = 0;
    {
        std::stringstream ss;
        iteration           = util::parseIteration( prims_path );
        size_t      it_loc  = prims_path.find("_it");
        std::string fname   = prims_path.substr( 0, it_loc );
        ss << fname << "_merged_it" << iteration << ".csv";
        o_path              = ss.str();
        io::savePrimitives<_PrimitiveT, typename PrimitiveMapT::mapped_type::const_iterator>( out_prims, o_path );
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
* \tparam _PointPrimitiveT     Wraps a point, exposing pos() and dir() functions. Concept: \ref rapter::PointPrimitive.
* \tparam _PrimitiveT          Wraps a primitive, exposing pos() and dir() functions. Concept: \ref rapter::LinePrimitive2.
* \tparam _PointContainerT     Holds the points. Concept: std::vector< \ref rapter::PointPrimitive >.
* \tparam _PrimitiveContainerT Holds the primitives grouped by GID in a two level structure. Concept: std::map< int, std::vector< \ref rapter::LinePrimitive2 > >.
* \tparam _Scalar              Floating point precision of primitives, points, etc. Concept: \ref rapter::PointPrimitive::Scalar.
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
                        , char                 const  mode
                        , int                         poplimit )
{
    //typedef rapter::MyPointPrimitiveDistanceFunctor _PointPrimitiveDistanceFunctor;
    typedef typename _PrimitiveContainerT::const_iterator   outer_const_iterator;
    //typedef typename outer_const_iterator::value_type::const_iterator inner_const_iterator;
    typedef typename _PointContainerT::value_type           PointT;
    typedef          Eigen::Matrix<_Scalar,3,1>             Position;

    typedef           std::vector< Position         >       ExtremaT;
    //typedef           std::map   < int, ExtremaT>              LidExtremaT;
    typedef           std::pair  < GidT   , LidT    >       GidLid;
    typedef           std::map   < GidLid, ExtremaT >       GidLidExtremaT;

    int err = EXIT_SUCCESS;

    bool changed = false;
    int haCount = 0, orphanReCount = 0;

    // Loop over all points, and select orphans
    do
    {
        changed = false;

        _PointPrimitiveDistanceFunctor distFunctor;

        // Populations
        GidPidVectorMap populations; // populations[gid] == std::vector<int> {pid0,pid1,...}
        if ( EXIT_SUCCESS == err )
        {
            err = processing::getPopulations( populations, points );
            CHECK( err, "getPopulations" );
        }

        // cache extrema
        GidLidExtremaT extremaMap;

        std::cout << "Orphan re-assigned ";
        for ( size_t pIdId = 0; pIdId != points.size(); ++pIdId )
        {
            //const PidT pId       = points[pIdId].getTag(_PointPrimitiveT::TAGS::PID);
            const GidT pointGId = points[pIdId].getTag(_PointPrimitiveT::TAGS::GID);
            typename _PrimitiveContainerT::const_iterator it = prims.find( pointGId );

            // try to find if at least one of the primitive of the group is big
            // if is not the case, we can potentially re-assign
            auto isBigPatch = [] (const _PrimitiveT& prim) { return prim.getTag(_PrimitiveT::TAGS::STATUS) != _PrimitiveT::STATUS_VALUES::SMALL; };

            if (    ( it == prims.end()                            )    // the patch this point is assigned to does not exist
                 || ( !containers::valueOf<_PrimitiveT>(it).size() )    // the patch this point is assigned to is empty (no primitives in it)
                 || ( std::find_if((*it).second.begin(), (*it).second.end(), isBigPatch) == (*it).second.end()) // there is no big patch in the group
               )  // the patch this point is assigned to is too small
            {

                // We here have an orphan, so we need to iterate over all primitives and get the closest distance < scale
                _Scalar  minDist = std::numeric_limits<_Scalar>::max();
                GidT     minGid  = -1;
                Position pos     = points[pIdId].pos();

                // now loop over all primitives to get the closest for
                for ( outer_const_iterator it1 = prims.begin(); it1 != prims.end(); ++it1 )
                {
                    LidT lid = 0; // primitive linear index in patch
                    for ( _inner_const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2, ++lid )
                    {
                        const GidT gid = (*it2).getTag(_PrimitiveT::TAGS::GID);

                        // if large patch found
                        if ( (*it2).getTag(_PrimitiveT::TAGS::STATUS) != _PrimitiveT::STATUS_VALUES::SMALL)
                        {
                            if ( extremaMap.find(GidLid(gid,lid)) == extremaMap.end() )
                            {
                                it2->template setExtentOutdated(); // we want to recalculate to be sure, points might have been reassigned
                                err = it2->template getExtent<_PointPrimitiveT>
                                        ( extremaMap[ GidLid(gid,lid) ]
                                        , points
                                        , scale
                                        , populations[gid].size() ? &(populations[gid]) : NULL );
                            }

                            _Scalar dist = distFunctor.eval( extremaMap[ GidLid(gid,lid) ], *it2, pos );

                            // store minimum distance
                            if ( dist < minDist )
                            {
                                if ( it2->getDistance(pos) < scale ) // added by Aron on 8/1/2015
                                {
                                    minDist = dist;
                                    minGid  = gid;
                                }
                                else
                                    ++haCount;
                            }
                        }
                    }
                }

                //std::cout << "Process orphan.... " << minDist << std::endl;
                if ( (minDist < scale) && (minDist >= _Scalar(0.)) )
                {
                    // reassign point
                    points[pIdId].setTag( _PointPrimitiveT::TAGS::GID, minGid );
                    ++orphanReCount;
                    if ( !(orphanReCount % 1000) )
                        std::cout << orphanReCount << " points, ";
                    //std::cout << pid << " " << minGid << ", ";
                    changed = true;
                }
            } //...if reassign point
        }//...for all points
    } while (changed);
    std::cout << std::endl;

    std::cout << "HA, did not add orphan to segment " << haCount << " times!" << std::endl;

    return EXIT_SUCCESS;
} //...adoptPoints()

/*! \brief Merges adjacent patches that have the same direction ID or are almost parallel.
 *
 *  \tparam _PatchPatchDistanceFunctorT  Concept: \ref rapter::RepresentativeSqrPatchPatchDistanceFunctorT.
 *  \param[in] patchPatchDistFunct       Distance functor between two patches, to define adjacency.
 *  \param[in] spatial_threshold         Two extrema should be at least this close to be merged. Concept: \ref MergeParams::spatial_threshold_mult == 3 * scale.
 *  \sa #decide_merge()
 */
template < class    _PrimitiveT
         , class    _PointPrimitiveT
         //, class    _inner_const_iterator
         , class    _InnerPrimitiveContainerT
         , class    _PrimitiveContainerT
         , class    _PointContainerT
         , typename _Scalar
         , class    _PatchPatchDistanceFunctorT
         , class    _PrimitiveDecideMergeFunctorT
         , class    _ComparedUidsT >
int Merging::mergeSameDirGids( _PrimitiveContainerT               & out_primitives
                             , _PointContainerT                   & points
                             , _PrimitiveContainerT        /*const&*/ primitives
                             , _Scalar                       const  scale
                             , _Scalar                       const  spatial_threshold
                             , _Scalar                       const  parallel_limit
                             , _PatchPatchDistanceFunctorT   const& /*patchPatchDistFunct*/
                             , _PrimitiveDecideMergeFunctorT const& primitiveDecideMergeFunct
                             , bool                                 preserveSmallPatches
                             , _ComparedUidsT                     & comparedUids
                             )
{
    typedef typename _PrimitiveContainerT::const_iterator      outer_const_iterator;
    typedef typename _InnerPrimitiveContainerT::const_iterator inner_const_iterator;
    typedef typename _PrimitiveContainerT::iterator            outer_iterator;
    typedef typename _InnerPrimitiveContainerT::iterator       inner_iterator;
    typedef           std::vector<Eigen::Matrix<_Scalar,3,1> > ExtremaT;
    typedef           std::map   < LidT, ExtremaT>             LidExtremaT;
    typedef           std::map   < GidT, LidExtremaT >         GidLidExtremaT;
    typedef           std::pair  < GidT, // map key
                                   LidT> // linear index in the array associated to the key
                      GidLid;
    typedef typename _ComparedUidsT::ElementT UidPairT;

    // Store the primitives that have been matched and must be ignored
    // First  (key)   = candidate,
    // Second (value) = reference.
    typedef std::set< LidT > IgnoreListT;
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
    UidT uid = 0; DidT maxDirGId = 0;
    GidLidExtremaT extrema; // <gid,lid> -> vector<x0, x1, ...>
    if ( EXIT_SUCCESS == err )
    {
        // for all patches
        for ( outer_iterator outer_it  = primitives.begin();
                                  (outer_it != primitives.end()); // we now handle error
                                 ++outer_it )
        {
            GidT gid = -2; // (-1 is "unset")
            LidT lid = 0;  // linear index of primitive in container (to keep track)
            // for all directions
            for ( inner_iterator inner_it  = containers::valueOf<_PrimitiveT>(outer_it).begin();
                                       (inner_it != containers::valueOf<_PrimitiveT>(outer_it).end());// we now handle error
                                      ++inner_it, ++lid )
            {
                //primitives.at(gid0).at(lid0).setTag(_PrimitiveT::USER_ID1, uid);
                inner_it->setTag( _PrimitiveT::USER_ID1, ++uid );
                maxDirGId = std::max( maxDirGId
                                    , static_cast<DidT>( inner_it->getTag(_PrimitiveT::TAGS::DIR_GID) )
                                    );

                // ignore small patches
                if ( inner_it->getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL )
                {
                    ignoreList.insert(gid);
                    continue;
                }

                // save patch gid
                if ( gid == -2 )
                {
                    gid = inner_it->getTag( _PrimitiveT::TAGS::GID );
                    // sanity check
                    if ( extrema.find(gid) != extrema.end() )   std::cerr << "[" << __func__ << "]: " << "GID not unique for patch...:-S" << std::endl;
                }

                if ( populations[gid].size() )
                {
                    err = inner_it->template getExtent<_PointPrimitiveT>
                                                               ( extrema[gid][lid]
                                                               , points
                                                               , scale
                                                               , &(populations[gid])
                                                               );
                }
                else
                    err = EXIT_FAILURE;

                if ( err != EXIT_SUCCESS )
                {
                    std::cerr << "Issue when computing extent of ("
                              << primitives.at(gid).at(lid).getTag(_PrimitiveT::TAGS::GID )     << ","
                              << primitives.at(gid).at(lid).getTag(_PrimitiveT::TAGS::DIR_GID ) << ")"
                              << std::endl
                              << "Ignored later... " << std::endl;

                    ignoreList.insert(gid);
                }
            } //...for primitives
        } //...for patches

        CHECK( err, "calcExtrema" )
    } //...getExtrema


    typedef typename GidLidExtremaT::const_iterator GidIt;
    typedef typename LidExtremaT::const_iterator    PrimIt;

#if 0
    // add a local uid attribute to primitives to be able to recognize them later
    DidT altMaxDirGid = 0;
    for ( GidIt gid_it = extrema.cbegin(); gid_it != extrema.cend(); ++gid_it )
    {
        const GidT gid0 = gid_it->first;
        for ( PrimIt prim_it = gid_it->second.cbegin(); prim_it != gid_it->second.cend(); ++prim_it, uid++ )
        {
            const GidT lid0 = std::distance<typename LidExtremaT::const_iterator>( gid_it ->second.cbegin(), prim_it );
            primitives.at(gid0).at(lid0).setTag(_PrimitiveT::USER_ID1, uid);
            altMaxDirGid = std::max( altMaxDirGid, static_cast<DidT>(primitives.at(gid0).at(lid0).getTag(_PrimitiveT::TAGS::DIR_GID)) );
        }
    }


    // estimate the maximum direction id, to be able to give new ids to refit primitives
    DidT max_dir_gid = 0;
    for ( outer_const_iterator outer_it  = primitives.begin(); (outer_it != primitives.end()); ++outer_it )
        for ( inner_const_iterator inner_it  = containers::valueOf<_PrimitiveT>(outer_it).begin();
                                   (inner_it != containers::valueOf<_PrimitiveT>(outer_it).end());
                                  ++inner_it )
        {
            max_dir_gid = std::max( max_dir_gid, static_cast<DidT>((*inner_it).getTag(_PrimitiveT::TAGS::DIR_GID)) );
        }
#endif
    std::cout << "[" << __func__ << "]: " << "max_dir_gid: " << maxDirGId << std::endl;

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

    //bool merged = false;

    // Reference traversal
    for ( GidIt gid_it = extrema.cbegin(); gid_it != extrema.cend(); ++gid_it )
    {
        const GidT gid0 = gid_it->first;

        for ( PrimIt prim_it = gid_it->second.cbegin(); prim_it != gid_it->second.cend(); ++prim_it )
        {
            // linear id of the reference
            const LidT lid0 = std::distance<typename LidExtremaT::const_iterator>( gid_it ->second.cbegin(), prim_it );

            // Reference primitive key, used to check if it has not been merged previously
            GidLid refGidLid (gid0, lid0);

            // check if this primitives has not been merged previously
            if (ignoreList.find(gid0) != ignoreList.end()) continue;

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

                GidT gid1 = gid_it1->first;

                // Given the ref. primitive, we check for merging candidates. They can be in:
                //   - the same map entry, after the ref. primitive in the linear Primitive array
                //   - the next map entries, any index
                for ( PrimIt prim_it1 = (gid_it1 == gid_it) ? ++PrimIt(prim_it) :       // for the same map entry, after the ref. primitive
                                                              gid_it1->second.cbegin(); // for the next map entries, any index
                      is0Valid && prim_it1 != gid_it1->second.cend();
                      ++prim_it1 )
                {

                    const LidT lid1 = std::distance<typename LidExtremaT::const_iterator>( gid_it1->second.cbegin(), prim_it1 );

                    // resume of the candidate primitive, used to record the matching operation for later and
                    // detect previous merge operations
                    GidLid candGidLid (gid1, lid1);

                    // Here we don't need to define
                    // bool is1Valid,
                    // calling continue is sufficient to jump to the next primitive after and merge,
                    // plus here check that a previous merge has not been recorded
                    if (ignoreList.find(gid1) != ignoreList.end()) continue;

                    const _PrimitiveT& prim1 = primitives.at(gid1).at(lid1);

                    PidT uid40 = prim0.getTag( _PrimitiveT::USER_TAGS::USER_ID4 ),
                         uid41 = prim1.getTag( _PrimitiveT::USER_TAGS::USER_ID4 );

                    UidPairT uid4Pair;
                    if ( uid40 > uid41 ) uid4Pair = UidPairT(uid41,uid40);
                    else                 uid4Pair = UidPairT(uid40,uid41);

                    if ( comparedUids.find( uid4Pair ) != comparedUids.end() )
                    {
                        comparedUids.incHits();
                        continue;
                    }


                    if (primitiveDecideMergeFunct.eval( prim_it->second,  // extrema 0
                                                        prim0,            // prim 0
                                                        prim_it1->second, // extrema 1
                                                        prim1,            // prim 1
                                                        scale))
                    {
                        //std::cout << " YES" << std::endl;

                        // record this to detect unmerged primitives later and invalidate both primitives
                        ignoreList.insert(gid0);
                        ignoreList.insert(gid1);

                        if (    ( prim0.getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL )
                             || ( prim1.getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL ) )
                        {
                            std::cout << "[" << __func__ << "]: " << "crap, small patches are merged..." << std::endl; fflush(stdout);
                            throw new std::runtime_error("asdf");
                        }

                        merging::merge( out_primitives,     // [out] Container storing merged primitives
                                        prim0,              // [in]  First primitive (can be invalidated during the call)
                                        populations[gid0],  // [in]  First primitive population (point ids)
                                        prim1,              // [in]  Second primitive
                                        populations[gid1],  // [in]  Second primitive population (point ids)
                                        points,             // [in]  Point cloud
                                        scale,              // [in]  Working scale (for refit)
                                        maxDirGId,          // [in,out] maximum direction id
                                        comparedUids.getMaxId() // [in,out] maximum new unique id
                                       );

                        //merged = true;

                        is0Valid = false;

                        //comparedUids.erase( uid4Pair );
                        comparedUids.eraseAny( uid4Pair );

                        continue;  // jump to the next couple ref/candidate primitive
                    }
                    else
                    {

                        auto uidPairIt = comparedUids.find( uid4Pair );
                        if ( uidPairIt != comparedUids.end() )
                        {
                            std::cout << "this shouldn't happen, why are we rechecking this pair: " << uid40 << "," << uid41 << std::endl;
                        }

                        comparedUids.insert( uid4Pair );
                    }
                    //else
                    //    std::cout << " NO" << std::endl;
                }
            }
        }
    }

    //typedef typename _PrimitiveContainerT::mapped_type::iterator inner_iterator;

    // we can now remove primitives that are not assigned to any points
    processing::eraseNonAssignedPrimitives<_PrimitiveT, inner_iterator>(out_primitives, points, preserveSmallPatches);

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
        if (l0.getTag(_PrimitiveT::TAGS::GID ) == l1.getTag(_PrimitiveT::TAGS::GID ))
            return l0.getTag(_PrimitiveT::TAGS::DIR_GID ) < l1.getTag(_PrimitiveT::TAGS::DIR_GID );
        return l0.getTag(_PrimitiveT::TAGS::GID ) < l1.getTag(_PrimitiveT::TAGS::GID );
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

            if (l0.getTag(_PrimitiveT::TAGS::DIR_GID ) != l1.getTag(_PrimitiveT::TAGS::DIR_GID ) ||
                l0.getTag(_PrimitiveT::TAGS::GID ) != l1.getTag(_PrimitiveT::TAGS::GID ) )
                return true;
        }
    }

    // if we reach that point, that means that input/output are similar
    return false;
} //...Merging::mergeSameDirGids()

} //...namespace rapter

#undef CHECK

#endif // RAPTER_MERGING_HPP
