#ifndef GF2_CANDIDATEGENERATOR_HPP
#define GF2_CANDIDATEGENERATOR_HPP

#include "globfit2/my_types.h"                      // PCLPointAllocator
#include "globfit2/util/util.hpp"                   // parseIteration()
#include "globfit2/processing/util.hpp"             // calcPopulations()
#include "globfit2/optimization/energyFunctors.h"   // MyPointPrimitiveDistanceFunctor
#include "globfit2/io/io.h"                         // readPrimities,savePrimitives,etc.
#include "globfit2/util/diskUtil.hpp"               // saveBackup
#include "globfit2/processing/angle_util.hpp"       // selectAngles

//debug
#include "pcl/point_types.h" // debug
#include "pcl/point_cloud.h" // debug

#define isSMALL(prim) (prim.getTag(_PrimitiveT::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL)
#define notSMALL(prim) (prim.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::SMALL)
#define isPROMOTED(gid,lid) (promoted.find( GidLid(gid,lid) ) != promoted.end())
#define notPROMOTED(gid,lid) (promoted.find( GidLid(gid,lid) ) == promoted.end())

namespace GF2
{
    namespace vis {
        //inline void
    }

    namespace cgen
    {
        template <typename _PrimitiveT, typename _PrimitiveContainerT>
        struct FilterChosenFunctor
        {
            typedef std::pair<int,int> GidLid;

            inline FilterChosenFunctor( _PrimitiveContainerT &out_prims, std::set<GidLid> const& chosen, std::set<GidLid> const& promoted )
                : _out_prims( out_prims )
                , _chosen( chosen )
                , _promoted( promoted )
            {
                if ( out_prims.size() )
                    std::cerr << "[" << __func__ << "]: " << "out_prims not empty...it was assumed to be empty..." << std::endl;
            }

            inline int eval( _PrimitiveT const& prim, int lid )
            {
                const int gid = prim.getTag(_PrimitiveT::GID);
                //const GidLid gidLid( gid, lid );

                // all promoted, and their derivatives fail here, since they have it set to the original lid
                bool keep   = prim.getTag(_PrimitiveT::USER_ID1) == _PrimitiveT::TAG_UNSET;
                bool demote = false; // flipped to true, if unchoosen promoted patch, that needs to be demoted to small, and put into output
                if ( !keep )
                {
                    // receiver needs to be chosen, if receiver was promoted
                    keep |= _chosen.find( GidLid(gid,prim.getTag(_PrimitiveT::USER_ID1)) ) != _chosen.end();

                    if ( !keep )
                    {
                        // the original patch was promoted, and not chosen, so we need to demote it, and put it back
//                        std::cout << "skip: gid:" << prim.getTag(_PrimitiveT::GID) << ", lid: " << lid
//                                  << ", tag1: " << prim.getTag(_PrimitiveT::USER_ID1)
//                                  << ", active: " << prim.getTag(_PrimitiveT::STATUS)
//                                  << std::endl;
                        if ( lid == prim.getTag(_PrimitiveT::USER_ID1) ) // this is the original promoted patch
                        {
                            demote = true;
                            keep = true;
                        }
                    }
//                    else
//                        std::cout << "KEEP: gid:" << prim.getTag(_PrimitiveT::GID) << ", lid: " << lid << ", tag1: " << prim.getTag(_PrimitiveT::USER_ID1) << std::endl;

                    if ( prim.getTag(_PrimitiveT::STATUS) != _PrimitiveT::STATUS_VALUES::UNSET )
                        std::cout << "[" << __func__ << "]: " << "status should not be " << prim.getTag(_PrimitiveT::STATUS) << ", it should be UNSET" << std::endl;
                }

                //bool patch_promoted = std::find_if( _promoted.begin(), _promoted.end(), [&gid](GidLid const& e){return e.first == gid;} ) != _promoted.end();
                //bool patch_chosen = std::find_if( _chosen.begin(), _chosen.end(), [&gid](GidLid const& e){return e.first == gid;} ) != _chosen.end();

                //if ( !patch_promoted || patch_chosen ) // add, if not promoted, or if promoted and chosen
                if ( keep )
                {
                    _PrimitiveT &added = containers::add( _out_prims, gid, prim );
                    if ( demote )
                        added.setTag( _PrimitiveT::STATUS, _PrimitiveT::STATUS_VALUES::SMALL );

                    return notSMALL( added );
                }

                return 0;
            }

            _PrimitiveContainerT &_out_prims;
            std::set<GidLid> const &_chosen, &_promoted;
        };
    }

    template <class _PrimitiveT, class _PrimitiveContainerT, class _GeneratedT, class _CopiedT, class _PromotedT>
    inline bool output( _PrimitiveT & cand, _PrimitiveContainerT &out_prims, _GeneratedT &generated, _CopiedT &copied, int &nlines
                      , _PromotedT const& promoted, int const gid, int const lid0, int const closest_angle_id, std::map<int,AnglesT> allowedAngles )
    {
        typedef std::pair<int,int> GidLid;
        typedef std::pair<int,int> DidAid;
#if 1
        // filter similar
        {
            for ( int i = 0; i != out_prims[gid].size(); ++i )
            {
                if ( cand.getTag(_PrimitiveT::DIR_GID) == out_prims[gid][i].getTag(_PrimitiveT::DIR_GID) )
                {
                    //float eps = (cand0.template dir() - out_prims[gid0][i].template dir()).array().abs().sum();
                    float ang = GF2::angleInRad( cand.template dir(), out_prims[gid][i].template dir() );
                    if ( ang < 1.e-6 )
                    {
                        //std::cout << "SIMILAR: " << cand0.toString() << " vs. " << out_prims[gid0][i].toString() << std::endl;
                        DidAid didAid( cand.getTag(_PrimitiveT::DIR_GID),closest_angle_id ); // direction id, closest angle id
                        copied[ cand.getTag(_PrimitiveT::GID) ].insert( didAid );

                        return false;
                    }
//                    else
//                        std::cout << "not similar: " << cand0.toString() << " vs. " << out_prims[gid0][i].toString() << std::endl;
                }
            }
        }
#endif

        // This is a new candidate, make sure formulate knows that
        cand.setTag( _PrimitiveT::STATUS, _PrimitiveT::STATUS_VALUES::UNSET );
        if ( isPROMOTED(gid,lid0) )
        {
            cand.setTag( _PrimitiveT::USER_ID1, lid0 ); // receiverLid

            // record responsibility
            ++generated[ GidLid(gid,lid0) ]; // receiver
            //++generated[ GidLid(gid1,lid1) ]; // sender
        }
        _PrimitiveT &added = containers::add( out_prims, gid, cand ); // insert into output
        ++nlines;                                  // keep track of output size

        // debug
        int tmp_size = copied[ cand.getTag(_PrimitiveT::GID) ].size();

        // keep track of instances
        DidAid didAid( cand.getTag(_PrimitiveT::DIR_GID),closest_angle_id ); // direction id, closest angle id
        copied[ cand.getTag(_PrimitiveT::GID) ].insert( didAid );

        // debug
        if ( copied[cand.getTag(_PrimitiveT::GID) ].size() == tmp_size )
            std::cerr << "[" << __func__ << "][" << __LINE__ << "]: NOOOOO insertion, should not happen"
                      << cand.getTag( _PrimitiveT::GID ) << ", " << cand.getTag( _PrimitiveT::DIR_GID ) << std::endl;

        return true;
    }



    /*! \brief Main functionality to generate lines from points.
     *
     *  \tparam _PointPrimitiveDistanceFunctorT Concept: \ref MyPointPrimitiveDistanceFunctor.
     *  \tparam _PrimitiveT                     Concept: \ref GF2::LinePrimitive2.
     *  \param[in] smallThresh  Decides, what primitive counts as small. Usually a multiple of scale (4x...0.1x)
     *  \param[in] var_limit    How many output variables we are allowing. Default: 0, meaning no limit.
     *  \post Produces primitives with STATUS tag UNSET(-1) (new primitives) or ACTIVE(2) (previously selected primitives)
     */
    template <  class       _PrimitivePrimitiveAngleFunctorT
              , class       _PointPrimitiveDistanceFunctorT
              , class       _PrimitiveT
              , class       _PrimitiveContainerT
              , class       _PointContainerT
              , typename    _Scalar> int
    CandidateGenerator::generate( _PrimitiveContainerT                   & out_prims
                                , _PrimitiveContainerT              /*const*/& in_lines //non-const, becuase some of the primitives get promoted to large patches
                                , _PointContainerT                  const& points   // non-const to be able to add group tags
                                , _Scalar                           const  scale
                                , std::vector<_Scalar>              const& angles
                                , CandidateGeneratorParams<_Scalar> const& params
                                , _Scalar                           const  smallThreshMult
                                , bool                              const  safe_mode_arg
                                , int                               const  var_limit )
    {
        // _________ typedefs _________

        typedef typename _PointContainerT::value_type                      _PointPrimitiveT;
        typedef typename _PrimitiveContainerT::const_iterator              outer_const_iterator;
        typedef typename _PrimitiveContainerT::iterator                    outer_iterator;
        typedef typename _PrimitiveContainerT::mapped_type                 InnerContainerT;
        //typedef typename outer_const_iterator::value_type::const_iterator inner_const_iterator;
        typedef typename InnerContainerT::const_iterator                   inner_const_iterator;
        typedef typename InnerContainerT::iterator                         inner_iterator;
        typedef          std::pair<int,int>                                GidLid;                  // uniquely identifies a primitive by first: gid, second: linear id in innerContainer (vector).
        typedef          std::pair<int,int>                                DidAid;
        typedef          std::map< GidLid, int >                           GeneratedMapT;
        typedef          std::pair<std::pair<int,int>,int>                 GeneratedEntryT;

        // cache, so that it can be turned off
        bool safe_mode = safe_mode_arg;
        bool debug = false;

        // _________ error checks _________

        if ( out_prims.size() ) std::cerr << "[" << __func__ << "]: " << "warning, out_lines not empty!" << std::endl;
        if ( params.patch_population_limit < 0 ) { std::cerr << "[" << __func__ << "]: " << "error, popfilter is necessary!!!" << std::endl; return EXIT_FAILURE; }
        if ( angles[0] != _Scalar(0.) )
            throw new CandidateGeneratorException("angles[0] = 0 always! we need it for parallel");
        if ( params.small_mode != CandidateGeneratorParams<_Scalar>::SmallPatchesMode::IGNORE )
            throw new CandidateGeneratorException("small option is out of order, use IGNORE, meaning ignore small");

        // _________ variables _________

        int                                 nlines = 0; // output size
        // filter already copied directions
        std::map< int, std::set<DidAid> >   copied; // [gid] = [ <dir0,angle_id0>, <dir0,angle_id1>, <dir2,angle_id0>, ... ]
        // modified angular similarity
        const _Scalar                       angle_limit( params.angle_limit / params.angle_limit_div );
        GeneratedMapT                       generated; // records, how many extra candidates the input primitive generated in the output


        // count patch populations
        GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        processing::getPopulations( populations, points );

        // debug
        typedef pcl::PointXYZRGB PCLPointT;
        typedef pcl::PointCloud<PCLPointT> PCLCloudT;
        typedef typename PCLCloudT::Ptr PCLCloudPtrT;
        PCLCloudPtrT cloud( new PCLCloudT() );
        _PointPrimitiveT::template toCloud<PCLCloudPtrT, _PointContainerT, PCLPointAllocator<_PointPrimitiveT::Dim> >
                ( cloud, points );

        // _________ (1) promotion _________

        // upgrade primitives to large, and copy large primitives to output
        // added on 21/09/2014 by Aron
        std::set<GidLid> promoted; // Contains primitives, that were small, but now are large (active)
        const _Scalar smallThresh = smallThreshMult * scale;
        {
            Eigen::Matrix<_Scalar,Eigen::Dynamic,1> spatialSignif(1,1); // cache variable

            // either all (patches), or none (second iteration) have to be unset
            int unset_count = 0, input_count = 0;
            // for all input primitives
            for ( outer_iterator outer_it0 = in_lines.begin(); outer_it0 != in_lines.end(); ++outer_it0 )
            {
                int gid = (*outer_it0).first;
                int lid = 0;
                for ( inner_iterator inner_it0 = (*outer_it0).second.begin(); inner_it0 != (*outer_it0).second.end(); ++inner_it0, ++lid )
                {
                    // cache vars
                    _PrimitiveT& prim         = *inner_it0;
                    int          prim_status  = prim      .getTag(_PrimitiveT::TAGS::STATUS);

                    // bookkeeping
                    if ( prim_status == _PrimitiveT::STATUS_VALUES::UNSET )
                        ++unset_count;

                    // if unset (first iteration) or small (second+ iteration ), threshold it
                    if (    (prim_status == _PrimitiveT::STATUS_VALUES::SMALL)      // promotable OR
                         || (prim_status == _PrimitiveT::STATUS_VALUES::UNSET) )    // first ranking has to be done, since this is first iteration patchify output
                    {
                        // Promote: check, if patch is large by now
                        if ( prim.getSpatialSignificance( spatialSignif, points, scale, &(populations[gid]))(0) >= smallThresh )
                        {
                            // store primitives, that have just been promoted to large from small
                            if ( (prim_status == _PrimitiveT::STATUS_VALUES::SMALL) )
                                promoted.insert( GidLid(gid,lid) );

                            prim.setTag( _PrimitiveT::STATUS, _PrimitiveT::STATUS_VALUES::UNSET ); // set to unset, so that formulate can distinguish between promoted and
                        }
                        // Demote, if first iteration
                        else if (prim_status == _PrimitiveT::STATUS_VALUES::UNSET) // this should only happen in the first iteration
                        {
                            //std::cout << "demoting" << spatialSignif(0) << std::endl;
                            prim.setTag( _PrimitiveT::STATUS, _PrimitiveT::STATUS_VALUES::SMALL );
                        }
                        else if (prim_status != _PrimitiveT::STATUS_VALUES::SMALL) // this should only happen in the first iteration
                        {
                            std::cout << "large patch( status: " << prim_status << "): " << spatialSignif(0) << std::endl;
                        }
//                        else
//                            std::cout << "not promoting size " << spatialSignif(0) << " < " << smallThresh << std::endl;

                        // update cache
                        prim_status = prim.getTag( _PrimitiveT::STATUS );
                    } //...if not large

                    // copy to output
                    {
                        // copy input - to keep CHOSEN tag entries, so that we can start second iteration from a selection
                        auto &added = containers::add( out_prims, gid, *inner_it0 );

                        if ( isPROMOTED(gid,lid) )
                        {
                            ++generated[ GidLid(gid,lid) ]; // receiver
                            added.setTag( _PrimitiveT::USER_ID1, lid ); // make sure, we can withdraw this later
                        }

                        // update output count
                        if ( notSMALL(prim) )
                            ++nlines;

                        // smalls are invisible in this round
                        if ( (prim_status == _PrimitiveT::STATUS_VALUES::ACTIVE) )
                        {
                            // store position-direction combination to avoid duplicates
                            copied[ gid ].insert( DidAid(inner_it0->getTag(_PrimitiveT::DIR_GID),0) );
                            //std::cout << "added to copied[" << gid << "]:" << inner_it0->getTag(_PrimitiveT::DIR_GID) << std::endl;
                        }
                    } //...copy to output

                    // bookkeeping
                    ++input_count;
                } //...for primitives in patch
            } //...for patches

            // check for mistakes - either all (patches), or none (second iteration) have to be unset
            if ( unset_count && (unset_count != input_count) )
                throw new CandidateGeneratorException("All tags have to be either unset (first iteration) or set to active/small (second iteration), cannot be partial here");

            std::cout << "[" << __func__ << "]: " << "promoted " << promoted.size() << " primitives" << std::endl;
            if ( !promoted.size() && safe_mode )
            {
                std::cout << "[" << __func__ << "]: " << "\n\nDISABLING safe mode, since no promoted\n" << std::endl;
                safe_mode = false;
            }
        } //...promote

        // _________ (2) statistics _________

        // 2.1 ANGLES
        AnglesT angle_gens_in_rad;
        deduceGenerators<_Scalar>( angle_gens_in_rad, angles );
        std::cout<<"angle_gens_in_rad:";for(size_t vi=0;vi!=angle_gens_in_rad.size();++vi)std::cout<<angle_gens_in_rad[vi]*180./M_PI<<" ";std::cout << "\n";
        // estimate direction cluster angles
        DirAngleMapT dirAngles;        // intermediate storage to collect dir angle distributions
        std::map<int,AnglesT> allowedAngles;      // final storage to store allowed angles
        {
            DirAnglesFunctorOuter<_PrimitiveT,_PrimitiveContainerT,AnglesT> outerFunctor( in_lines, angles, /* output reference */ dirAngles );
            processing::filterPrimitives<_PrimitiveT,typename _PrimitiveContainerT::mapped_type>( in_lines, outerFunctor );

            if ( !dirAngles.size() )
                std::cerr << "[" << __func__ << "]: dirAngles.size(): " << dirAngles.size() << std::endl;

            selectAngles( allowedAngles, dirAngles, angles, angle_gens_in_rad );
        }

        // 2.2 Neighbourhoods
        //_Scalar dist = MyFinitePrimitiveToFinitePrimitiveCompatFunctor::eval( ex1, p1, ex2, p2 );

        // Status of primitives at this stage are:
        //  ACTIVE      for large primitives to use
        //  TAG_UNSET   for promoted primitives to give directions to
        //  SMALL       for primitives to ignore

        //debug
        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );

        // _________ (4) generation _________
        int gid0, gid1, dir_gid0, dir_gid1, lid0, lid1;
        gid0 = gid1 = dir_gid0 = dir_gid1 = _PrimitiveT::TAG_UNSET; // group tags cached
        for ( outer_const_iterator outer_it0  = in_lines.begin(); outer_it0 != in_lines.end(); ++outer_it0 )
        {
            gid0 = -2; // -1 is unset, -2 is unread
            lid0 =  0;
            for ( inner_const_iterator inner_it0  = (*outer_it0).second.begin(); inner_it0 != (*outer_it0).second.end(); ++inner_it0, ++lid0 )
            {
                // cache outer primitive
                _PrimitiveT const& prim0 = *inner_it0;
                dir_gid0                 = prim0.getTag( _PrimitiveT::DIR_GID );

                // cache group id of patch at first member
                     if ( gid0 == -2               )  gid0 = prim0.getTag( _PrimitiveT::GID ); // store gid of first member in patch
                else if ( gid0 != outer_it0->first )  std::cerr << "[" << __func__ << "]: " << "Not good, prims under one gid don't have same GID..." << std::endl;

                for ( outer_const_iterator outer_it1  = outer_it0; outer_it1 != in_lines.end(); ++outer_it1 )
                {
                    gid1 = -2; // -1 is unset, -2 is unread
                    lid1 = 0;
                    for ( inner_const_iterator inner_it1  = (outer_it0 == outer_it1) ? ++inner_const_iterator( inner_it0 )
                                                                                     : (*outer_it1).second.begin();
                                               inner_it1 != (*outer_it1).second.end();
                                             ++inner_it1, ++lid1 )
                    {
                        // cache inner primitive
                        _PrimitiveT const& prim1 = containers::valueOf<_PrimitiveT>( inner_it1 );
                        dir_gid1                 = prim1.getTag( _PrimitiveT::DIR_GID );

                        // cache group id of patch at first member
                             if ( gid1 == -2               )    gid1 = prim1.getTag( _PrimitiveT::GID );
                        else if ( gid1 != outer_it1->first )    std::cerr << "[" << __func__ << "]: " << "Not good, prims under one gid don't have same GID in inner loop..." << std::endl;

                        bool add0 = false, // add0: new primitive at location of prim0, with direction from prim1.
                             add1 = false; // add1: new primitive at location of prim1, with direction from prim0

#                       warning "This uses bias for all patches, fix later by calling this method with receive_all in the last iteration and fixing the contents there"

                        // prim0 (pos) needs to be active or promoted to receive directions
                        // prim1 (dir) needs to be active and not promoted to send direction
                        add0 = notSMALL(prim0) && notSMALL(prim1) && notPROMOTED(gid1,lid1);
                        // changed by Aron 08:15, 26/09/2014 to prevent many candidates, trust originals, don't add anything to them, just add to promoted ones
                        if ( safe_mode )    add0 &= isPROMOTED(gid0,lid0); // add cand0 only, if prim0 was promoted

                        // prim0 (dir) needs to be active and not promoted to send direction
                        // prim1 (pos) needs to be active or promoted to receive directions
                        add1 = notSMALL(prim0) && notSMALL(prim1) && notPROMOTED(gid0,lid0);
                        // changed by Aron 08:15, 26/09/2014 to prevent many candidates, trust originals, don't add anything to them, just add to promoted ones
                        if ( safe_mode )    add1 &= isPROMOTED(gid1,lid1); // add cand1 only, if prim1 was promoted

                        // find best rotation id and value
                        int closest_angle_id0 = 0, closest_angle_id1 = 0;
                        if ( add0 || add1 ) // speedup by skipping angle calculation, if won't add anyway... // added by Aron 12:46 18 Oct 2014
                        {
                            bool isAllowedAngle0 = allowedAngles.find(dir_gid1) != allowedAngles.end(), // do I have restrictions, when I copy dir1 to pos0 with this angle?
                                 isAllowedAngle1 = allowedAngles.find(dir_gid0) != allowedAngles.end(); // do I have restrictions, when I copy dir0 to pos1 with this angle?

                            // ANG0
                            // angdiff0 determines, if dir1 from pos1 can be copied to pos0 using allowed angles for dir1
                            _Scalar angdiff0 = _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>( prim0, prim1,
                                                                                                 isAllowedAngle0 ? allowedAngles[dir_gid1] : angles,
                                                                                                 &closest_angle_id0 );
                            // retarget closest angle id0 to angles instead of allowedAngles[dir_gid1]
                            if ( isAllowedAngle0 )
                                for ( int aa = 0; aa != angles.size(); ++aa )
                                    if ( angles[aa] == allowedAngles[dir_gid1][closest_angle_id0] )
                                    {
                                        closest_angle_id0 = aa;
                                        break;
                                    }

                            // ANG1
                            // angdiff1 determines, if dir0 from pos0 can be copied to pos1 using allowed angles for dir0
                            _Scalar angdiff1 = _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>( prim0, prim1,
                                                                                                 isAllowedAngle1 ? allowedAngles[dir_gid0] : angles,
                                                                                                 &closest_angle_id1 );
                            // retarget closest angle id1 to angles instead of allowedAngles[dir_gid0]
                            if ( isAllowedAngle1 )
                                for ( int aa = 0; aa != angles.size(); ++aa )
                                    if ( angles[aa] == allowedAngles[dir_gid0][closest_angle_id1] )
                                    {
                                        closest_angle_id1 = aa;
                                        break;
                                    }

                            // decide based on rotation difference (RECEIVE_SIMILAR)
                            //bool    close_ang  = angdiff < angle_limit;
                            add0              &= angdiff0 < angle_limit; // prim0 needs to be close to prim1 in angle
                            add1              &= angdiff1 < angle_limit; // prim1 needs to be close to prim0 in angle
                        }

                        // store already copied pairs
                        {
                            if ( copied[gid0].find(DidAid(dir_gid1,closest_angle_id0)) != copied[gid0].end() )
                                add0 = false;

                            if ( copied[gid1].find(DidAid(dir_gid0,closest_angle_id1)) != copied[gid1].end() )
                                add1 = false;
                        }

                        if ( !add0 && !add1 )
                            continue;

                        // debug
                        if (    ((dir_gid0 == 127) || (dir_gid1 == 127))
                             && (     (((gid0==71) || (gid0==80) || (gid0==64) || (gid0==67)) && add0)
                                  ||  (((gid1==71) || (gid1==80) || (gid1==64) || (gid1==67)) && add1)
                                )
                           )
                        {
                            std::cout << "sotp" << std::endl;
                            debug = true;
                        }

                        bool added0 = add0, added1 = add1;

                        AnglesT single_gen0, single_gen1;
                        genAngles( single_gen0, angles[closest_angle_id0], angle_gens_in_rad );
                        genAngles( single_gen1, angles[closest_angle_id1], angle_gens_in_rad );

                        // copy line from pid to pid1
                        _PrimitiveT cand0;
                        if ( add0 && prim0.generateFrom(cand0, prim1, closest_angle_id0, angles, _Scalar(1.)) )
                        {
                            added0 = output( cand0, out_prims, generated, copied, nlines
                                           , promoted, gid0, lid0, closest_angle_id0, allowedAngles );

                            if ( allowedAngles.find(dir_gid0) == allowedAngles.end() && single_gen0.size() )
                            {
                                processing::appendAnglesFromGenerators( allowedAngles[dir_gid0], single_gen0, /* no_paral: */ false, true, true );
                            }
                        } //...add0

                        _PrimitiveT cand1;
                        if ( add1 && prim1.generateFrom(cand1, prim0, closest_angle_id1, angles, _Scalar(-1.)) ) // changed to 1 on 22 12 2014
                        {
                            added1 = output( cand1, out_prims, generated, copied, nlines
                                           , promoted, gid1, lid1, closest_angle_id1, allowedAngles );

                            if ( allowedAngles.find(dir_gid1) == allowedAngles.end() && single_gen1.size() )
                            {
                                processing::appendAnglesFromGenerators( allowedAngles[dir_gid1], single_gen1, /* no_paral: */ false, true, true );
                            }
                        } //...if add1

                        //debug
                        {
                            if ( debug )
                            {
                                vptr->setBackgroundColor( .9, .9, .9 );
                                vptr->addPointCloud( cloud );

                                char line_name[255];
                                sprintf( line_name, "prim0" );
                                _PrimitiveT::template draw<_PointPrimitiveT>( /*   primitive: */ prim0
                                                                              , /*      points: */ points
                                                                              , /*   threshold: */ scale
                                                                              , /*     indices: */ &(populations[gid0])
                                                                              , /*      viewer: */ vptr
                                                                              , /*   unique_id: */ line_name
                                                                              , /*      colour: */ 1,0,0
                                                                              , /* viewport_id: */ 0
                                                                              );
                                sprintf( line_name, "prim1" );
                                _PrimitiveT::template draw<_PointPrimitiveT>( /*   primitive: */ prim1
                                                                              , /*      points: */ points
                                                                              , /*   threshold: */ scale
                                                                              , /*     indices: */ &(populations[gid1])
                                                                              , /*      viewer: */ vptr
                                                                              , /*   unique_id: */ line_name
                                                                              , /*      colour: */ 0,1,0
                                                                              , /* viewport_id: */ 0
                                                                              );
                                if ( added0 )
                                {
                                    sprintf( line_name, "cand0" );
                                    _PrimitiveT::template draw<_PointPrimitiveT>( /*   primitive: */ cand0
                                                                                  , /*      points: */ points
                                                                                  , /*   threshold: */ scale
                                                                                  , /*     indices: */ &(populations[gid0])
                                                                                  , /*      viewer: */ vptr
                                                                                  , /*   unique_id: */ line_name
                                                                                  , /*      colour: */ 0,0,1
                                                                                  , /* viewport_id: */ 0
                                                                                  );
                                }
                                if ( added1 )
                                {
                                    sprintf( line_name, "cand1" );
                                    _PrimitiveT::template draw<_PointPrimitiveT>( /*   primitive: */ cand1
                                                                                  , /*      points: */ points
                                                                                  , /*   threshold: */ scale
                                                                                  , /*     indices: */ &(populations[gid1])
                                                                                  , /*      viewer: */ vptr
                                                                                  , /*   unique_id: */ line_name
                                                                                  , /*      colour: */ 0,1,1
                                                                                  , /* viewport_id: */ 0
                                                                                  );
                                }
                                std::cout << "gid: " << gid0 << ", dir_gid: " << dir_gid0
                                          << ", gid1: " << gid1 << ", dir_gid1: " << dir_gid1
                                          << std::endl;

                                if ( added0 )
                                {
                                    std::cout << "copied[" << gid0 << "]: ";
                                    for ( auto it = copied[gid0].begin(); it != copied[gid0].end(); ++it )
                                        std::cout << "<" << it->first << "," << it->second << "," << angles[it->second] * 180./M_PI << ">, ";
                                    std::cout << "\twas looking for " << "<" << dir_gid1 << "," << closest_angle_id0 << "," << angles[closest_angle_id0] * 180./M_PI << ">";
                                    std::cout << std::endl;
                                }

                                if ( added1 )
                                {
                                    std::cout << "copied[" << gid1 << "]: ";
                                    for ( auto it = copied[gid1].begin(); it != copied[gid1].end(); ++it )
                                        std::cout << "<" << it->first << "," << it->second << "," << angles[it->second] * 180./M_PI << ">, ";
                                    std::cout << "\twas looking for " << "<" << dir_gid0 << "," << closest_angle_id1 << "," << angles[closest_angle_id1] * 180./M_PI << ">";
                                    std::cout << std::endl;
                                }

                                if ( debug && (added0 || added1) )
                                {
                                    vptr->spin();
                                    debug = false;
                                }
                                vptr->removeAllShapes();
                            }
                            else
                                std::cout << "gid: " << gid0 << ", dir_gid: " << dir_gid0 << std::endl;
                        } //...debug

                    } //...for l3
                } //...for l2
            } //...for l1
        } //...for l0

        // /home/bontius/workspace/globOpt/data/scenes/triangle_60d$ ../glob_opt --generate -sc 0.04 -al 0.4 -ald 1 --small-mode 0 --patch-pop-limit 25 --angle-gens 60 -p patches.csv --assoc points_primitives.csv --small-thresh-mult 1
        // new: /home/bontius/workspace/globOpt/data/scenes/pearl1_test$ ../glob_opt --generate -sc 0.02 -al 0.6 -ald 1 --small-mode 0 --patch-pop-limit 20 --angle-gens 90 -p primitives_merged_it0.csv --assoc points_primitives_it0.csv --small-thresh-mult 0.05

        // ___________ LIMIT VARIABLES _______________
        int ret = EXIT_SUCCESS;
        if ( (var_limit > 0) && (nlines > var_limit) )
        {
            int active_count = nlines; // records how many variables actives generated among themselves, these we cannot filter

            // (1) sort promoted patches descending based on how many variables they generated
            // in: generated, containing < <gid,lid> , generated count> entries for promoted locations only (no actives)
            // out: ranks, containing a sorted list of generated
            std::vector< GeneratedEntryT > ranks;
            for ( auto it = generated.begin(); it != generated.end(); ++it )
            {
                // in_lines[it.first.first][it.first.second].getTag( _PrimitiveT::STATUS ) )
                //if ( promoted.find( GidLid(it->first) ) != promoted.end() )

                ranks.push_back( *it );
                active_count -= it->second;

                // was promoted
                // std::cout << "[" << __func__ << "]: " << "promoted[" << it->first.first << "," << it->first.second << "] generated " << it->second << std::endl;
            }

            // don't do this, if the actives alone were exceeding the limit...
            if ( active_count >= var_limit )
            {
                std::cerr << "[" << __func__ << "]: " << "!!!!!!!! ALL " << nlines << " ACTIVE or ACTIVEs " << active_count << " > " << var_limit << "var_limit, cannot limit vars! !!!!!!!" << std::endl;
            }
            else
            {
                // sort descending:
                std::sort( ranks.begin(), ranks.end(), [](GeneratedMapT::value_type const& a, GeneratedMapT::value_type const& b){ return a.second > b.second;} );

                // (2) Select entries, that get to stay promoted
                // in: ranks, containing a descending list of variable counts
                // out: chosen, a list of promoted locations, that should stay promoted in out_prims, others should be demoted back
                std::cout << "var_limit: " << var_limit << ", active: " << active_count << ", nlines: " << nlines << std::endl;
                std::set< GidLid > chosen;
                {
                    int vars = active_count;
                    for ( int i = 0; i != ranks.size(); ++i )
                    {
                        vars += ranks[i].second;
                        if ( vars <= var_limit )
                        {
                            chosen.insert( ranks[i].first );
                        }
                    }
                }

                // (3) Filter out_prims based on chosen
                //
                auto tmp = out_prims; out_prims.clear();
                cgen::FilterChosenFunctor<_PrimitiveT, _PrimitiveContainerT> functor( out_prims, chosen, promoted );
                nlines = processing::filterPrimitives<_PrimitiveT, InnerContainerT>( tmp, functor );
                // return remaining variable count
                ret = ranks.size() - chosen.size();
            } // if too many actives already
        } // filter

        // log
        std::cout << "[" << __func__ << "]: " << "finished generating, we now have " << nlines << " candidates" << std::endl;

        return ret;
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
            generatorParams.safe_mode = pcl::console::find_switch( argc, argv, "--safe-mode" );
            if ( generatorParams.safe_mode )
                std::cout << "[" << __func__ << "]: " << "__________________________\n__________________________RUNNING SAFE______________________\n_____________________________" << std::endl;
            pcl::console::parse_argument( argc, argv, "--var-limit", generatorParams.var_limit );
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
                if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
                {
                    std::cout << "[" << __func__ << "]: " << "Usage:\t " << argv[0] << " --generate \n";
                    std::cout << "\t --cloud " << cloud_path << "\n";
                    std::cout << "\t -sc,--scale " << generatorParams.scale << "\n";
                    std::cout << "\t -p,--prims" << input_prims_path << "\n";
                    std::cout << "\t -a,--assoc" << associations_path << "\n";
                    std::cout << "\t --small-thresh-mult " << generatorParams.small_thresh_mult << "\t get's multiplied to scale, and serves as big/small threshold for primitives\n";

                    // linkage mode (full_min, full_max, squared_min, repr_min)
                    std::cout << "\t [--mode *" << generatorParams.printPatchDistMode() << "*\t";
                    for ( size_t m = 0; m != mode_opts.size(); ++m )
                        std::cout << "|" << mode_opts[m];
                    std::cout << "]\n";

                    std::cout << "\t [-al,--angle-limit " << generatorParams.angle_limit << "]\n";
                    std::cout << "\t [-ald,--angle-limit-div " << generatorParams.angle_limit_div << "]\n";
                    std::cout << "\t [--patch-dist-limit " << generatorParams.patch_dist_limit_mult << "]\n";
                    std::cout << "\t [--angle-gens "; for(size_t vi=0;vi!=angle_gens.size();++vi)std::cout<<angle_gens[vi]<<","; std::cout << "]\n";
                    std::cout << "\t [--patch-pop-limit " << generatorParams.patch_population_limit << "]\n";
                    std::cout << "\t [--small-mode " << generatorParams.small_mode << "\t | 0: IGNORE, 1: RECEIVE_SIMILAR, 2: RECEIVE_ALL]\n";
                    std::cout << "\t [--no-paral]\n";
                    std::cout << "\t [--safe-mode]\n";
                    std::cout << "\t [--var-limit " << generatorParams.var_limit << "\t Decides how many variables we want as output. 0 means unlimited.]\n";
                    std::cout << std::endl;

                    return EXIT_FAILURE;
                }
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

        bool no_paral = pcl::console::find_switch( argc, argv, "--no-paral");
        // Read desired angles
        if ( EXIT_SUCCESS == err )
        {
            processing::appendAnglesFromGenerators( generatorParams.angles, angle_gens, no_paral, true );
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
        int ret = EXIT_SUCCESS;
        if ( EXIT_SUCCESS == err )
        {
            ret = CandidateGenerator::generate< MyPrimitivePrimitiveAngleFunctor, MyPointPrimitiveDistanceFunctor, _PrimitiveT >
                                              ( primitives, patches, points
                                                , generatorParams.scale
                                                , generatorParams.angles
                                                , generatorParams
                                                , generatorParams.small_thresh_mult
                                                , generatorParams.safe_mode
                                                , generatorParams.var_limit );

            //if ( err != EXIT_SUCCESS ) std::cerr << "[" << __func__ << "]: " << "generate exited with error! Code: " << err << std::endl;
            if ( ret != 0 )
                std::cout << "not all patches were promoted( " << ret << " left), will need rerun on same threshold..." << std::endl;
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

        return std::max(ret,err);
    } // ...CandidateGenerator::generateCli()

} // ... ns GF2

#undef isSMALL
#undef notSMALL
#undef isPROMOTED
#undef notPROMOTED

#endif // GF2_CANDIDATEGENERATOR_HPP
