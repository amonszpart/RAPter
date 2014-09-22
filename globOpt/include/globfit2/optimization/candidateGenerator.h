#ifndef __GF2_CANDIDATEGENERATOR_H__
#define __GF2_CANDIDATEGENERATOR_H__

#include "globfit2/parameters.h"                         // CandidateGeneratorParams

namespace GF2
{
    class CandidateGenerator
    {
        public:
             /*! \brief                  Step 1. Generates primitives from a cloud. Reads "cloud.ply" and saves "candidates.csv".
              *  \param argc             Contains --cloud cloud.ply, and --scale scale.
              *  \param argv             Contains --cloud cloud.ply, and --scale scale.
              *  \return                 EXIT_SUCCESS.
              */
            template < class    _PrimitiveContainerT
                     , class    _PointContainerT
                     , typename _Scalar
                     , class    _PointPrimitiveT
                     , class    _PrimitiveT
                     >
            static inline int
            generateCli( int argc, char** argv );

            /*! \brief Main functionality to generate lines from points.
             *
             *  \tparam _PointPrimitiveDistanceFunctorT Concept: \ref MyPointPrimitiveDistanceFunctor.
             *  \tparam _PrimitiveT                     Concept: \ref GF2::LinePrimitive2.
             *  \param[in] smallThresh Decides, what primitive counts as small. Usually a multiple of scale (4x...0.1x)
             *  \post Produces primitives with STATUS tag UNSET(-1) (new primitives) or ACTIVE(2) (previously selected primitives)
             */
            template <  class       PrimitivePrimitiveAngleFunctorT // concept: energyFunctors.h::PrimitivePrimitiveAngleFunctor
                      , class       _PointPrimitiveDistanceFunctorT
                      , class       _PrimitiveT
                      , class       PrimitiveContainerT             // concept: std::vector<std::vector<LinePrimitive2>>
                      , class       PointContainerT                 // concept: std::vector<PointPrimitive>
                      , typename    Scalar
                      >
            static inline int
            generate( PrimitiveContainerT                   &  out_lines
                    , PrimitiveContainerT                   &  in_lines // non-const, because some of the primitives are promoted
                    , PointContainerT                  const&  points
                    , Scalar                           const   scale
                    , std::vector<Scalar>              const&  angles
                    , CandidateGeneratorParams<Scalar> const&  params
                    , Scalar                           const  smallThresh );

    }; //...class CandidateGenerator
} // ...ns::GF2

//_____________________________________________________________________________________________________________________
// HPP
#include "globfit2/my_types.h"                      // PCLPointAllocator
#include "globfit2/util/util.hpp"                   // parseIteration()
#include "globfit2/processing/util.hpp"             // calcPopulations()
#include "globfit2/optimization/energyFunctors.h"   // MyPointPrimitiveDistanceFunctor
#include "globfit2/io/io.h"                         // readPrimities,savePrimitives,etc.
#include "globfit2/util/diskUtil.hpp"               // saveBackup

namespace GF2
{
    inline bool equal2D( int l0, int l1, int l2, int l3 ) { return (l0 == l2) && (l1 == l3); }

    /*! \brief Called from \ref CandidateGenerator::generate() via \ref processing::filterPrimitives, decides if a primitive has the GID looked for.
     * \tparam _PrimitiveT Concept: \ref GF2::LinePrimtive2.
     */
    template <class _PointPrimitiveDistanceFunctor, class _PrimitiveT, class _PointPrimitiveT, typename _Scalar>
    struct NearbyPrimitivesFunctor
    {
            NearbyPrimitivesFunctor( _PointPrimitiveT const& point, _Scalar scale ) : _point(point), _scale(scale) {}

            inline int eval( _PrimitiveT const& prim, int const lid )
            {
                // store if inside scale
                if ( _PointPrimitiveDistanceFunctor::template eval<_Scalar>(_point, prim) < _scale )
                {
                    _gidLids.push_back( std::pair<int,int>(prim.getTag(_PrimitiveT::GID), lid) );
                } //...if dist

                return 0;
            }

            _PointPrimitiveT                  _point;
            _Scalar                           _scale;
            std::vector< std::pair<int,int> > _gidLids; //!< Output unique IDs of primitives.
    }; //NearbyPrimitivesFunctor


    /*! \brief Called from getOrphanGids via \ref processing::filterPrimitives, decides if a primitive has the GID looked for.
     */
    template <class _PrimitiveT>
    struct FindGidFunctor
    {
            FindGidFunctor( int gid ) : _gid(gid) {}
            inline int eval( _PrimitiveT const& prim, int /*lid*/ ) const { return prim.getTag(_PrimitiveT::GID) == _gid; }
            int  _gid; //!< \brief GID looked for.
    }; //...FindGidFunctor

    /*! \brief Gathers groupIDs that don't have primitives selected for them.
     *
     *  \tparam _PrimitiveT           Concept: \ref GF2::LinePrimtive2.
     *  \tparam _inner_const_iterator Concept: _PrimitiveContainerT::value_type::const_iterator if vector, _PrimitiveContainerT::mapped_type::const_iterator if map.
     *  \tparam _PidContainerT        Container holding the orphan point ids. Concept: std::set<int>.
     *  \tparam _GidContainerT        Container holding the patch ids (GIDs) with orphan points. Concept: std::set<int>.
     *  \tparam _PointContainerT      Concept: vector<_PointPrimitiveT>.
     *  \tparam _PrimitiveContainerT  Concept: vector< vector< _PrimitiveT> >
     *
     * \param[in] point_gid_tag      Identifies the field, where the GID is stored in the pointPrimitives. Concept: _PointPrimitiveT::GID
     */
    template < class _PrimitiveT
             , class _inner_const_iterator
             , class _PidContainerT
             , class _GidContainerT
             , class _PointContainerT
             , class _PrimitiveContainerT
             >
    inline int getOrphanGids( _GidContainerT            & orphan_gids
                            , _PointContainerT     const& points
                            , _PrimitiveContainerT const& prims
                            , int                  const  point_gid_tag
                            , _PidContainerT            * orphan_pids   = NULL )
    {
        // select unassigned points
        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            int gid = points[pid].getTag( point_gid_tag );

            FindGidFunctor<_PrimitiveT> functor(gid);
            if ( !processing::filterPrimitives< _PrimitiveT, _inner_const_iterator>
                                              ( prims, functor )                                )
            {
                orphan_gids.insert( gid );
                if ( orphan_pids )
                    orphan_pids->insert( pid );
            }
        }
        return EXIT_SUCCESS;
    }


    template <  class       _PrimitivePrimitiveAngleFunctorT
              , class       _PointPrimitiveDistanceFunctorT
              , class       _PrimitiveT
              , class       _PrimitiveContainerT
              , class       _PointContainerT
              , typename    _Scalar> int
    CandidateGenerator::generate( _PrimitiveContainerT                   & out_lines
                                , _PrimitiveContainerT              /*const*/& in_lines //non-const, becuase some of the primitives get promoted to large patches
                                , _PointContainerT                  const& points   // non-const to be able to add group tags
                                , _Scalar                           const  scale
                                , std::vector<_Scalar>              const& angles
                                , CandidateGeneratorParams<_Scalar> const& params
                                , _Scalar                           const  smallThreshMult )
    {
        typedef typename _PointContainerT::value_type                      _PointPrimitiveT;
        typedef typename _PrimitiveContainerT::const_iterator              outer_const_iterator;
        typedef typename _PrimitiveContainerT::iterator                    outer_iterator;
        //typedef typename outer_const_iterator::value_type::const_iterator inner_const_iterator;
        typedef typename _PrimitiveContainerT::mapped_type::const_iterator inner_const_iterator;
        typedef typename _PrimitiveContainerT::mapped_type::iterator       inner_iterator;
        typedef          std::pair<int,int>                                GidLid;                  // uniquely identifies a primitive by first: gid, second: linear id in innerContainer (vector).

        if ( out_lines.size() ) std::cerr << "[" << __func__ << "]: " << "warning, out_lines not empty!" << std::endl;
        if ( params.patch_population_limit < 0 ) { std::cerr << "[" << __func__ << "]: " << "error, popfilter is necessary!!!" << std::endl; return EXIT_FAILURE; }

        int nlines = 0; // output size
        // filter already copied directions
        std::map< int, std::set<int> > copied; // [gid] = [ dir0, dir 1, dir2, ... ]

        // (2) Mix and Filter
        // count patch populations
        GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        {
            processing::getPopulations( populations, points );
        }

        // upgrade primitives to large,
        // and copy large primitives to output
        // added on 21/09/2014 by Aron
        std::set<GidLid> promoted; // Contains primitives, that were small, but now are large (active)
        const _Scalar smallThresh = smallThreshMult * scale;
        {
            Eigen::Matrix<_Scalar,Eigen::Dynamic,1> spatialSignif(1,1); // cache variable

            // either all (patches), or none (second iteration) have to be unset
            int unset_count = 0, input_count = 0;
            // for all input primitives
            for ( outer_iterator outer_it0  = in_lines.begin(); outer_it0 != in_lines.end(); ++outer_it0 )
            {
                int gid = (*outer_it0).first;
                int lid = 0;
                for ( inner_iterator inner_it0  = (*outer_it0).second.begin(); inner_it0 != (*outer_it0).second.end(); ++inner_it0, ++lid )
                {
                    // cache vars
                    _PrimitiveT& prim         = *inner_it0;
                    int          prim_status  = prim.getTag(_PrimitiveT::TAGS::STATUS);

                    // bookkeeping
                    if ( prim_status == _PrimitiveT::STATUS_VALUES::UNSET )
                        ++unset_count;

                    // if unset (first iteration) or small (second+ iteration ), threshold it
                    if (    (prim_status == _PrimitiveT::STATUS_VALUES::SMALL)      // promotable OR
                         || (prim_status == _PrimitiveT::STATUS_VALUES::UNSET) )    // first ranking has to be done, since this is first iteration patchify output
                    {
                        // Promote: check, if patch is large by now
                        if ( prim.getSpatialSignificance( spatialSignif, points, scale)(0) > smallThresh )
                        {
                            std::cout << "promoting " << spatialSignif(0) << std::endl;
                            // store primitives, that have just been promoted to large from small
                            if ( (prim_status == _PrimitiveT::STATUS_VALUES::SMALL) )
                                promoted.insert( GidLid(gid,lid) );

                            prim.setTag( _PrimitiveT::STATUS, _PrimitiveT::STATUS_VALUES::UNSET ); // set to unset, so that formulate can distinguish between promoted and
                        }
                        // Demote, if first iteration
                        else if (prim_status == _PrimitiveT::STATUS_VALUES::UNSET) // this should only happen in the first iteration
                        {
                            std::cout << "demoting" << spatialSignif(0) << std::endl;
                            prim.setTag( _PrimitiveT::STATUS, _PrimitiveT::STATUS_VALUES::SMALL );

                        }
                        else
                        {
                            std::cout << "large patch( status: " << prim_status << "): " << spatialSignif(0) << std::endl;
                        }

                        // update cache
                        prim_status = prim.getTag( _PrimitiveT::STATUS );
                    } //...if not large

#warning This does not copy promoted, remember to account for them later
                    // copy to output
                    //  if first iteration: copy UNSET
                    //  if second iteration: copy ACTIVE and SMALL, there should be no unset, except promoted, which we don't want to copy
                    //if ( (prim_status == _PrimitiveT::STATUS_VALUES::ACTIVE) || (prim_status == _PrimitiveT::STATUS_VALUES::SMALL) ) // don't copy promoted...
                    {
                        // copy input - to keep CHOSEN tag entries, so that we can start second iteration from a selection
                        containers::add( out_lines, inner_it0->getTag( _PrimitiveT::GID ), *inner_it0 );

                        // smalls are invisible in this round
                        if ( (prim_status == _PrimitiveT::STATUS_VALUES::ACTIVE) )
                        {
                            // store position-direction combination to avoid duplicates
                            copied[ inner_it0->getTag(_PrimitiveT::GID) ].insert( inner_it0->getTag(_PrimitiveT::DIR_GID) );
                            // update output count
                            ++nlines;
                        }
                    } //...if active, copy to output

                    // bookkeeping
                    ++input_count;
                } //...for primitives in patch
            } //...for patches

            // check for mistakes - either all (patches), or none (second iteration) have to be unset
            if ( unset_count && (unset_count != input_count) )
                throw new std::runtime_error("All tags have to be either unset (first iteration) or set to active/small (second iteration), cannot be partial here");
//            if ( unset_count == input_count ) // first iteration
//                promoted.clear();

            std::cout << "[" << __func__ << "]: " << "promoted " << promoted.size() << " primitives" << std::endl;
        } //...sort size

        // convert local fits
        const _Scalar angle_limit( params.angle_limit / params.angle_limit_div );

        // The following cycles go through each patch and each primitive in each patch and pair them up with each other patch and the primitives within.
        // So, 00 - 00, 00 - 01, 00 - 02, ..., 01 - 01, 01 - 02, 01 - 03, ... , 10 - 10, 10 - 11, 10 - 12, ..., 20 - 20, 20 - 21,  ... (ab - cd), etc.
        // for a = 0 {   for b = 0 {   for c = a {   for d = b {   ...   } } } }

        // depr:
#if 0
        // copy input to output (to keep chosen tags)
        for ( outer_const_iterator outer_it0  = in_lines.begin(); outer_it0 != in_lines.end(); ++outer_it0 )
        {
            for ( inner_const_iterator inner_it0  = (*outer_it0).second.begin(); inner_it0 != (*outer_it0).second.end(); ++inner_it0 )
            {
                // copy input - to keep status tag entries, so that we can start second iteration from a selection
                containers::add( out_lines, inner_it0->getTag( _PrimitiveT::GID ), *inner_it0 );
                // store position-direction combination to avoid duplicates
                copied[ inner_it0->getTag(_PrimitiveT::GID) ].insert( inner_it0->getTag(_PrimitiveT::DIR_GID) );
                // update output count
                ++nlines;
            }
        }
#endif

        // Status of primitives at this stage are:
        //  ACTIVE      for large primitives to use
        //  TAG_UNSET   for promoted primitives to give directions to
        //  SMALL       for primitives to ignore

        int gid0, gid1, dir_gid0, dir_gid1, lid0, lid1;
        gid0 = gid1 = dir_gid0 = dir_gid1 = _PrimitiveT::TAG_UNSET; // group tags cached
        // OUTER0 (a)
        for ( outer_const_iterator outer_it0  = in_lines.begin(); outer_it0 != in_lines.end(); ++outer_it0 )
        {
            // INNER1 (b)
            gid0 = -2; // -1 is unset, -2 is unread
            lid0 = 0;
            for ( inner_const_iterator inner_it0  = (*outer_it0).second.begin(); inner_it0 != (*outer_it0).second.end(); ++inner_it0, ++lid0 )
            {
                // cache outer primitive
                _PrimitiveT const& prim0 = *inner_it0;
                dir_gid0                 = prim0.getTag( _PrimitiveT::DIR_GID );

                // cache group id of patch at first member
                if ( gid0 == -2 )
                    gid0 = prim0.getTag( _PrimitiveT::GID ); // store gid of first member in patch
                else if ( (*outer_it0).first != gid0 )
                    std::cerr << "[" << __func__ << "]: " << "Not good, prims under one gid don't have same GID..." << std::endl;

                // OUTER1 (c)
                for ( outer_const_iterator outer_it1  = outer_it0; outer_it1 != in_lines.end(); ++outer_it1 )
                {
                    gid1 = -2; // -1 is unset, -2 is unread
                    lid1 = 0;
                    // INNER1 (d)
                    for ( inner_const_iterator inner_it1  = (outer_it0 == outer_it1) ? ++inner_const_iterator( inner_it0 )
                                                                                     : (*outer_it1).second.begin();
                                               inner_it1 != (*outer_it1).second.end();
                                             ++inner_it1, ++lid1 )
                    {
                        // cache inner primitive
                        _PrimitiveT const& prim1 = containers::valueOf<_PrimitiveT>( inner_it1 );
                        dir_gid1                 = prim1.getTag( _PrimitiveT::DIR_GID );

                        // cache group id of patch at first member
                        if ( gid1 == -2 )
                            gid1 = prim1.getTag( _PrimitiveT::GID );
                        else if ( (*outer_it1).first != gid1 )
                            std::cerr << "[" << __func__ << "]: " << "Not good, prims under one gid don't have same GID in inner loop..." << std::endl;

                        bool add0 = false, // add0: new primitive at location of prim0, with direction from prim1.
                             add1 = false; // add1: new primitive at location of prim1, with direction from prim0

                        // find best rotation id and value
                        int     closest_angle_id = 0;
                        //_Scalar closest_angle    = _Scalar( 0 );
                        {
                            _Scalar angdiff    = _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>( prim0, prim1, angles, &closest_angle_id );
                            //closest_angle      = angles[ closest_angle_id ];

                            // decide based on rotation difference (RECEIVE_SIMILAR)
                            bool    close_ang  = angdiff < angle_limit;
                            add0              |= close_ang; // prim0 needs to be close to prim1 in angle
                            add1              |= close_ang; // prim1 needs to be close to prim0 in angle
                        }

                        if ( params.small_mode == CandidateGeneratorParams<_Scalar>::SmallPatchesMode::RECEIVE_ALL ) // unused branch
                        {
                            add0 |=    ( prim0.getTag(_PrimitiveT::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL  )
                                    && ( prim1.getTag(_PrimitiveT::STATUS) == _PrimitiveT::STATUS_VALUES::ACTIVE ); // prim0 needs to be small to copy the dir of prim1.
                            // |= (populations[gid0].size() < params.patch_population_limit);
                            add1 |=    ( prim1.getTag(_PrimitiveT::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL )
                                    && ( prim0.getTag(_PrimitiveT::STATUS) == _PrimitiveT::STATUS_VALUES::ACTIVE ); // prim1 needs to be small to copy the dir of prim0.
                            // |= (populations[gid1].size() < params.patch_population_limit)
                            std::cerr << "[" << __func__ << "]: " << "SMALL_MODE==RECEIVE_ALL____________________ARE YOU SURE______________________???" << std::endl;
                        }
                        else if ( params.small_mode == CandidateGeneratorParams<_Scalar>::SmallPatchesMode::IGNORE )
                        {
#warning "This uses bias for small patches, fix later by calling this method with receive_all in the last iteration and fixing the contents there"

                            // cand0 = pos0, dir1
                            add0 &= (   (prim0.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::SMALL)       // prim0 (pos) needs to be active or promoted to receive directions
                                     && (prim1.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::SMALL)       // prim1 (dir) needs to be active             to send    directions
                                     && promoted.find( GidLid(gid1,lid1) ) == promoted.end() // prim1 needs to be active and not promotoed to send direction
                                    );

                            add1 &= (   (prim0.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::SMALL )      // prim0 (dir) needs to be active             to send    directions
                                     && (prim1.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::SMALL ) // prim1 (pos) needs to be active or promoted to receive directions
                                     && promoted.find( GidLid(gid0,lid0) ) == promoted.end() // prim0 needs to be active and not promotoed to send direction
                                    );
                        }

                        // store already copied pairs
                        {
                            if ( copied[gid0].find(dir_gid1) != copied[gid0].end() )
                                add0 = false;

                            if ( copied[gid1].find(dir_gid0) != copied[gid1].end() )
                                add1 = false;
                        }

                        if ( !add0 && !add1 )
                            continue;

                        // copy line from pid to pid1
                        _PrimitiveT cand0;
                        if ( add0 && prim0.generateFrom(cand0, prim1, closest_angle_id, angles, _Scalar(-1.)) )
                        {
                            // prepare
                            //_PrimitiveT cand0 = _PrimitiveT( prim0.pos(), dir0 );
                            //cand0.setTag( _PrimitiveT::GID    , prim0.getTag(_PrimitiveT::GID) );
                            //cand0.setTag( _PrimitiveT::DIR_GID, prim1.getTag(_PrimitiveT::DIR_GID) ); // recently changed this from GID

                            // debug - show here, that promote worked
                            {
                                if ( (prim0.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::ACTIVE) || (prim1.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::ACTIVE) )
                                    std::cout << "possible promote. cand0.pos = prim0(" << prim0.getTag(_PrimitiveT::STATUS) << ").pos, cand0.dir = prim1(" << prim1.getTag(_PrimitiveT::STATUS)<< ").dir()\n";
                            }

                            // This is a new candidate, make sure formulate knows that
                            cand0.setTag( _PrimitiveT::STATUS, _PrimitiveT::STATUS_VALUES::UNSET );

                            // insert
                            int check = containers::add( out_lines, gid0, cand0 ).getTag( _PrimitiveT::GID );
                            ++nlines;

                            // debug
                            if ( check != cand0.getTag(_PrimitiveT::GID) ) std::cerr << "gid tag copy check failed" << std::endl;
                            int tmp_size = copied[ cand0.getTag(_PrimitiveT::GID) ].size();

                            // keep track of instances
                            copied[ cand0.getTag(_PrimitiveT::GID) ].insert( cand0.getTag(_PrimitiveT::DIR_GID) );

                            // debug
                            if ( copied[cand0.getTag(_PrimitiveT::GID)].size() == tmp_size )
                                std::cerr << "[" << __func__ << "][" << __LINE__ << "]: NOOOOO insertion, should not happen"
                                          << cand0.getTag( _PrimitiveT::GID ) << ", " << cand0.getTag( _PrimitiveT::DIR_GID )
                                          << std::endl;
                        }

                        _PrimitiveT cand1;
                        if ( add1 && prim1.generateFrom(cand1, prim0, closest_angle_id, angles, _Scalar(1.)) )
                        {
                            //_PrimitiveT cand1 = _PrimitiveT( prim1.pos(), dir1 );
                            //_PrimitiveT cand1( prim1, prim0, -closest_angle );
                            //cand1.setTag( _PrimitiveT::GID    , prim1.getTag(_PrimitiveT::GID    ) );
                            //cand1.setTag( _PrimitiveT::DIR_GID, prim0.getTag(_PrimitiveT::DIR_GID) );

                            // debug - show here, that promote worked
                            {
                                if ( (prim0.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::ACTIVE) || (prim1.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::ACTIVE) )
                                    std::cout << "possible promote. cand1.pos = prim1(" << prim1.getTag(_PrimitiveT::STATUS) << ").pos, cand1.dir = prim0(" << prim0.getTag(_PrimitiveT::STATUS)<< ").dir()\n";
                            }

                            // This is a new candidate, make sure formulate knows that
                            cand1.setTag( _PrimitiveT::STATUS, _PrimitiveT::STATUS_VALUES::UNSET );

                            // copy line from pid1 to pid
                            int check = containers::add( out_lines, gid1, cand1 ).getTag( _PrimitiveT::DIR_GID );
                            ++nlines;

                            // debug
                            if ( check != cand1.getTag(_PrimitiveT::DIR_GID) ) std::cerr << "dirgid tag copy check failed" << std::endl;
                            int tmp_size = copied[ cand1.getTag(_PrimitiveT::GID) ].size();

                            // keep track of instances
                            copied[ cand1.getTag(_PrimitiveT::GID) ].insert(  cand1.getTag(_PrimitiveT::DIR_GID) );

                            // debug
                            if ( copied[cand1.getTag(_PrimitiveT::GID)].size() == tmp_size )
                                std::cerr << "[" << __func__ << "][" << __LINE__ << "]: NOOOOO insertion, should not happen for " << cand1.getTag( _PrimitiveT::GID ) << ", " << cand1.getTag( _PrimitiveT::DIR_GID ) << std::endl;
                        } //...if add1
                    } //...for l3
                } //...for l2
            } //...for l1
        } //...for l0

        // log
        std::cout << "[" << __func__ << "]: " << "finished generating, we now have " << nlines << " candidates" << std::endl;

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::generate()

    //! \brief                  Step 1. Generates primitives from a cloud. Reads "cloud.ply" and saves "candidates.csv".
    //! \param argc             Contains --cloud cloud.ply, and --scale scale.
    //! \param argv             Contains --cloud cloud.ply, and --scale scale.
    //! \return                 EXIT_SUCCESS.
    template < class    _PrimitiveContainerT
             , class    _PointContainerT
             , typename _Scalar
             , class    _PointPrimitiveT
             , class    _PrimitiveT
             >
    int
    CandidateGenerator::generateCli( int    argc
                                   , char** argv )
    {
        typedef typename _PrimitiveContainerT::value_type InnerPrimitiveContainerT;
        //typedef typename PointContainerT::value_type PointPrimitiveT;
        int err = EXIT_SUCCESS;

        CandidateGeneratorParams<Scalar> generatorParams;
        std::string                 cloud_path              = "./cloud.ply";
        std::vector<Scalar>         angle_gens              = { Scalar(90.) };
        std::string                 mode_string             = "representative_sqr";
        std::vector<std::string>    mode_opts               = { "representative_sqr" };
        std::string                 input_prims_path        = "patches.csv";
        std::string                 associations_path       = "points_primitives.csv";

        // parse input
        if ( err == EXIT_SUCCESS )
        {
            bool valid_input = true;

            // cloud
            if ( (pcl::console::parse_argument( argc, argv, "--cloud", cloud_path) < 0)
                 && !boost::filesystem::exists( cloud_path ) )
            {
                std::cerr << "[" << __func__ << "]: " << "--cloud does not exist: " << cloud_path << std::endl;
                valid_input = false;
            }

            // scale
            if ( (pcl::console::parse_argument( argc, argv, "--scale", generatorParams.scale) < 0) && (pcl::console::parse_argument( argc, argv, "-sc", generatorParams.scale) < 0) )
            {
                std::cerr << "[" << __func__ << "]: " << "--scale is compulsory" << std::endl;
                valid_input = false;
            }

            if (    (pcl::console::parse_argument( argc, argv, "-p", input_prims_path) < 0)
                 && (pcl::console::parse_argument( argc, argv, "--prims", input_prims_path) < 0)
                 && (!boost::filesystem::exists(input_prims_path)) )
            {
                std::cerr << "[" << __func__ << "]: " << "-p or --prims is compulsory" << std::endl;
                valid_input = false;
            }

            if ( (pcl::console::parse_argument( argc, argv, "--small-thresh-mult", generatorParams.small_thresh_mult) < 0) )
            {
                std::cerr << "[" << __func__ << "]: " << "--small-thresh-mult is compulsory" << std::endl;
                valid_input = false;
            }

            if (    (pcl::console::parse_argument( argc, argv, "-a", associations_path) < 0)
                 && (pcl::console::parse_argument( argc, argv, "--assoc", associations_path) < 0)
                 && (!boost::filesystem::exists(associations_path)) )
            {
                std::cerr << "[" << __func__ << "]: " << "-a or --assoc is compulsory" << std::endl;
                valid_input = false;
            }

            pcl::console::parse_argument( argc, argv, "--angle-limit", generatorParams.angle_limit );
            pcl::console::parse_argument( argc, argv, "-al", generatorParams.angle_limit );
            pcl::console::parse_argument( argc, argv, "--angle-limit-div", generatorParams.angle_limit_div );
            pcl::console::parse_argument( argc, argv, "-ald", generatorParams.angle_limit_div );
            pcl::console::parse_argument( argc, argv, "--patch-dist-limit", generatorParams.patch_dist_limit_mult ); // gets multiplied by scale
            pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens );
            pcl::console::parse_argument( argc, argv, "--patch-pop-limit", generatorParams.patch_population_limit );

            // patchDistMode
            pcl::console::parse_argument( argc, argv, "--mode", mode_string );
            generatorParams.parsePatchDistMode( mode_string );
            // refit
            if ( pcl::console::find_switch( argc, argv, "--patch-refit" ) )
            {
                std::cerr << "[" << __func__ << "]: " << "--patch-refit option has been DEPRECATED. exiting." << std::endl;
                return EXIT_FAILURE;
            }
            //pcl::console::parse_argument( argc, argv, "--patch-refit", patch_refit_mode_string );
            //generatorParams.parseRefitMode( patch_refit_mode_string );

            // small_mode
            {
                int small_mode = 0;
                pcl::console::parse_argument( argc, argv, "--small-mode", small_mode );
                generatorParams.small_mode = static_cast<CandidateGeneratorParams<Scalar>::SmallPatchesMode>( small_mode );
            }

            // print usage
            {
                std::cerr << "[" << __func__ << "]: " << "Usage:\t " << argv[0] << " --generate \n";
                std::cerr << "\t --cloud " << cloud_path << "\n";
                std::cerr << "\t -sc,--scale " << generatorParams.scale << "\n";
                std::cerr << "\t -p,--prims" << input_prims_path << "\n";
                std::cerr << "\t -a,--assoc" << associations_path << "\n";
                std::cerr << "\t --small-thresh-mult " << generatorParams.small_thresh_mult << "\t get's multiplied to scale, and serves as big/small threshold for primitives\n";

                // linkage mode (full_min, full_max, squared_min, repr_min)
                std::cerr << "\t [--mode *" << generatorParams.printPatchDistMode() << "*\t";
                for ( size_t m = 0; m != mode_opts.size(); ++m )
                    std::cerr << "|" << mode_opts[m];
                std::cerr << "]\n";

                std::cerr << "\t [-al,--angle-limit " << generatorParams.angle_limit << "]\n";
                std::cerr << "\t [-ald,--angle-limit-div " << generatorParams.angle_limit_div << "]\n";
                std::cerr << "\t [--patch-dist-limit " << generatorParams.patch_dist_limit_mult << "]\n";
                std::cerr << "\t [--angle-gens "; for(size_t vi=0;vi!=angle_gens.size();++vi)std::cerr<<angle_gens[vi]<<","; std::cerr << "]\n";
                std::cerr << "\t [--patch-pop-limit " << generatorParams.patch_population_limit << "]\n";
                std::cerr << "\t [--small-mode " << generatorParams.small_mode << "\t | 0: IGNORE, 1: RECEIVE_SIMILAR, 2: RECEIVE_ALL]\n";
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
            processing::appendAnglesFromGenerators( generatorParams.angles, angle_gens, true );
        } //...read angles

        // Read points
        PointContainerT points;
        if ( EXIT_SUCCESS == err )
        {
            err = io::readPoints<_PointPrimitiveT>( points, cloud_path );
            if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
        } //...read points

        std::vector<std::pair<int,int> > points_primitives;
        io::readAssociations( points_primitives, associations_path, NULL );
        for ( size_t i = 0; i != points.size(); ++i )
        {
            // store association in point
            points[i].setTag( PointPrimitiveT::GID, points_primitives[i].first );
        }

        // read primitives
        typedef std::map<int, InnerPrimitiveContainerT> PrimitiveMapT;
        _PrimitiveContainerT initial_primitives;
        PrimitiveMapT patches;
        {
            std::cout << "[" << __func__ << "]: " << "reading primitives from " << input_prims_path << "...";
            io::readPrimitives<_PrimitiveT, InnerPrimitiveContainerT>( initial_primitives, input_prims_path, &patches );
            std::cout << "reading primitives ok (#: " << initial_primitives.size() << ")\n";
        } //...read primitives

        //_____________________WORK_______________________
        //_______________________________________________

        // Generate
        //PrimitiveContainerT primitives;
        PrimitiveMapT primitives;
        if ( EXIT_SUCCESS == err )
        {
            err = CandidateGenerator::generate< MyPrimitivePrimitiveAngleFunctor, MyPointPrimitiveDistanceFunctor, _PrimitiveT >
                                              ( primitives, patches, points, generatorParams.scale, generatorParams.angles, generatorParams
                                                , generatorParams.small_thresh_mult);

            if ( err != EXIT_SUCCESS ) std::cerr << "[" << __func__ << "]: " << "generate exited with error! Code: " << err << std::endl;
        } //...generate

    #if 0 // these shouldn't change here
        // Save point GID tags
        if ( EXIT_SUCCESS == err )
        {
            std::string assoc_path = boost::filesystem::path( cloud_path ).parent_path().string() + "/" + "points_primitives.csv";

            util::saveBackup( assoc_path );
            err = io::writeAssociations<PointPrimitiveT>( points, assoc_path );

            if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "saveBackup or writeAssociations exited with error! Code: " << err << std::endl;
            else                        std::cout << "[" << __func__ << "]: " << "wrote to " << assoc_path << std::endl;

        } //...save Associations
    #endif

        // save primitives
        std::string o_path = boost::filesystem::path( cloud_path ).parent_path().string() + "/";
        if ( EXIT_SUCCESS == err )
        {
            std::string output_prims_path( o_path + "candidates.csv" );
            {
                int iteration = 0;
                iteration = util::parseIteration( input_prims_path ) + 1;
                std::stringstream ss;
                ss << o_path << "candidates_it" << iteration << ".csv";
                output_prims_path = ss.str();
            }

            util::saveBackup( output_prims_path );
            err = io::savePrimitives<_PrimitiveT,typename InnerPrimitiveContainerT::const_iterator>( /* what: */ primitives, /* where_to: */ output_prims_path );

            if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "saveBackup or savePrimitive exited with error! Code: " << err << std::endl;
            else                        std::cout << "[" << __func__ << "]: " << "wrote to " << output_prims_path << std::endl;
        } //...save primitives

        return err;
    } // ...CandidateGenerator::generateCli()

} // ... ns GF2


#endif // __GF2_CANDIDATEGENERATOR_H__

