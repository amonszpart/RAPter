#ifndef GF2_CANDIDATEGENERATOR_HPP
#define GF2_CANDIDATEGENERATOR_HPP

#include <map>

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
#include "pcl/console/parse.h"

#define isSMALL(prim) (prim.getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL)
#define notSMALL(prim) (prim.getTag(_PrimitiveT::TAGS::STATUS) != _PrimitiveT::STATUS_VALUES::SMALL)
#define isPROMOTED(gid,lid) (promoted.find( GidLid(gid,lid) ) != promoted.end())
#define notPROMOTED(gid,lid) (promoted.find( GidLid(gid,lid) ) == promoted.end())

namespace GF2
{
    namespace cgen
    {
        template <typename _PrimitiveT, typename _PrimitiveContainerT>
        struct FilterChosenFunctor
        {
            typedef std::pair<GidT,LidT> GidLid;

            inline FilterChosenFunctor( _PrimitiveContainerT &out_prims, std::set<GidLid> const& chosen, std::set<GidLid> const& promoted )
                : _out_prims( out_prims )
                , _chosen( chosen )
                , _promoted( promoted )
            {
                if ( out_prims.size() )
                    std::cerr << "[" << __func__ << "]: " << "out_prims not empty...it was assumed to be empty..." << std::endl;
            }

            inline LidT eval( _PrimitiveT const& prim, LidT lid )
            {
                const GidT gid = prim.getTag(_PrimitiveT::TAGS::GID);
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
//                        std::cout << "skip: gid:" << prim.getTag(_PrimitiveT::TAGS::GID) << ", lid: " << lid
//                                  << ", tag1: " << prim.getTag(_PrimitiveT::USER_ID1)
//                                  << ", active: " << prim.getTag(_PrimitiveT::TAGS::STATUS)
//                                  << std::endl;
                        if ( lid == prim.getTag(_PrimitiveT::USER_ID1) ) // this is the original promoted patch
                        {
                            demote = true;
                            keep = true;
                        }
                    }
//                    else
//                        std::cout << "KEEP: gid:" << prim.getTag(_PrimitiveT::TAGS::GID) << ", lid: " << lid << ", tag1: " << prim.getTag(_PrimitiveT::USER_ID1) << std::endl;

                    if ( prim.getTag(_PrimitiveT::TAGS::STATUS) != _PrimitiveT::STATUS_VALUES::UNSET )
                        std::cout << "[" << __func__ << "]: " << "status should not be " << prim.getTag(_PrimitiveT::TAGS::STATUS) << ", it should be UNSET" << std::endl;
                }

                //bool patch_promoted = std::find_if( _promoted.begin(), _promoted.end(), [&gid](GidLid const& e){return e.first == gid;} ) != _promoted.end();
                //bool patch_chosen = std::find_if( _chosen.begin(), _chosen.end(), [&gid](GidLid const& e){return e.first == gid;} ) != _chosen.end();

                //if ( !patch_promoted || patch_chosen ) // add, if not promoted, or if promoted and chosen
                if ( keep )
                {
                    _PrimitiveT &added = containers::add( _out_prims, gid, prim );
                    if ( demote )
                        added.setTag( _PrimitiveT::TAGS::STATUS, _PrimitiveT::STATUS_VALUES::SMALL );

                    return notSMALL( added );
                }

                return 0;
            }

            _PrimitiveContainerT &_out_prims;
            std::set<GidLid> const &_chosen, &_promoted;
        };
    }

    template <class _PrimitiveT, class _PrimitiveContainerT, class _GeneratedT, class _CopiedT, class _PromotedT>
    inline bool output( _PrimitiveT & cand, _PrimitiveContainerT &out_prims, _GeneratedT &generated, _CopiedT &copied, LidT &nlines
                      , _PromotedT const& promoted, GidT const gid, LidT const lid0, int const closest_angle_id, typename _PrimitiveT::Scalar const genAngle )
    {
        typedef typename _PrimitiveT::Scalar Scalar;

        // filter similar
        for ( LidT i = 0; i != out_prims[gid].size(); ++i )
        {
            if ( cand.getTag(_PrimitiveT::TAGS::DIR_GID) == out_prims[gid][i].getTag(_PrimitiveT::TAGS::DIR_GID) )
            {
                //float eps = (cand0.template dir() - out_prims[gid0][i].template dir()).array().abs().sum();
                Scalar ang = GF2::angleInRad( cand.template dir(), out_prims[gid][i].template dir() );
                if ( ang < 1.e-6 )
                {
                    //std::cout << "SIMILAR: " << cand0.toString() << " vs. " << out_prims[gid0][i].toString() << std::endl;
                    DidAid didAid( cand.getTag(_PrimitiveT::TAGS::DIR_GID),closest_angle_id ); // direction id, closest angle id
                    copied[ cand.getTag(_PrimitiveT::TAGS::GID) ].insert( didAid );

                    return false;
                }
                //                    else
                //                        std::cout << "not similar: " << cand0.toString() << " vs. " << out_prims[gid0][i].toString() << std::endl;
            } //...if same dId
        } //...for all prims in this patch (gId)

        // This is a new candidate, make sure formulate knows that
        cand.setTag( _PrimitiveT::TAGS::STATUS, _PrimitiveT::STATUS_VALUES::UNSET );
        // Make sure variable scheduling can find this candidate based on lid0
        if ( isPROMOTED(gid,lid0) )
        {
            cand.setTag( _PrimitiveT::USER_ID1, lid0 ); // receiverLid

            // record responsibility
            ++generated[ GidLid(gid,lid0) ]; // receiver
            //++generated[ GidLid(gid1,lid1) ]; // sender
        }
        // GEN_ANGLE
        if ( genAngle != _PrimitiveT::GEN_ANGLE_VALUES::UNSET )
            cand.setTag( _PrimitiveT::TAGS::GEN_ANGLE, genAngle );
        // Add
        /*_PrimitiveT &added = */containers::add( out_prims, gid, cand ); // insert into output
        ++nlines;                                                     // keep track of output size

        // debug
        LidT tmp_size = copied[ cand.getTag(_PrimitiveT::TAGS::GID) ].size();

        // keep track of instances
        DidAid didAid( cand.getTag(_PrimitiveT::TAGS::DIR_GID),closest_angle_id ); // direction id, closest angle id
        copied[ cand.getTag(_PrimitiveT::TAGS::GID) ].insert( didAid );

        // debug
        if ( copied[cand.getTag(_PrimitiveT::TAGS::GID) ].size() == tmp_size )
            std::cerr << "[" << __func__ << "][" << __LINE__ << "]: NOOOOO insertion, should not happen"
                      << cand.getTag( _PrimitiveT::TAGS::GID ) << ", " << cand.getTag( _PrimitiveT::TAGS::DIR_GID ) << std::endl;

        return true;
    }

    /*! \brief Adds a candidate at location prim0, with direction from prim1.
     *  \param[in/out] allowedAngles
     *  \param[in/out] copied           [gid] = [ <dir0,angle_id0>, <dir0,angle_id1>, <dir2,angle_id0>, ... ]
     *  \param[in/out] generated        Records, how many extra candidates the input primitive generated in the output.
     *  \param[in/out] nLines           Keeps track of overall output size.
     *  \param[in/out] aliases          Keeps track of diretion ids and their assigned generator angles. If not set, not enabled. NULL is used at the second phase, when aliases are added.
     */
    template < class _PrimitivePrimitiveAngleFunctorT, class _AliasesT
             , class _PrimitiveT, typename _Scalar, class _AnglesT, class _PromotedT
             , class _AllowedAnglesT, class _CopiedT, class _GeneratedT, class _PrimitiveContainerT, class _PointContainerT>
    inline int addCandidate( _PrimitiveT        const& prim0
                           , _PrimitiveT        const& prim1
                           , LidT               const  lid0
                           , LidT               const  lid1
                           , int                const  safe_mode
                           , int                const  allow_promoted
                           , _Scalar            const  angle_limit
                           , _AnglesT           const& angles
                           , _AnglesT           const& angle_gens_in_rad
                           , _PromotedT         const& promoted
                           , _AllowedAnglesT         & allowedAngles
                           , _CopiedT                & copied
                           , _GeneratedT             & generated
                           , LidT                    & nLines
                           , _PrimitiveContainerT    & out_prims
                           , _PointContainerT   const& points
                           , Scalar             const  scale
                           , _AliasesT               * aliases
                           , bool               const  tripletSafe  = false
                           , bool               const  verbose      = false
                           )
    {
        typedef typename _PointContainerT::value_type PointPrimitiveT;

        const GidT gid0     = prim0.getTag( _PrimitiveT::TAGS::GID );
        const GidT gid1     = prim1.getTag( _PrimitiveT::TAGS::GID );
        const DidT dir_gid0 = prim0.getTag( _PrimitiveT::TAGS::DIR_GID );
        const DidT dir_gid1 = prim1.getTag( _PrimitiveT::TAGS::DIR_GID );

        if ( gid0 == 50 || gid1 == 50 )
            std::cout << "at " << gid0 << "," << dir_gid0 << " from " << gid1  << ", " << dir_gid1 << std::endl;

        bool add0 = false;

        add0 = notSMALL(prim0) && notSMALL(prim1);
        if ( !allow_promoted )
            add0 &= notPROMOTED(gid1,lid1); // sender cannot be promoted
        // changed by Aron 08:15, 26/09/2014 to prevent many candidates, trust originals, don't add anything to them, just add to promoted ones
        if ( safe_mode )
            add0 &= isPROMOTED(gid0,lid0); // add cand0 only, if prim0 was promoted

        if ( !add0 )
            return false;

        int closest_angle_id0 = 0;
        // do I have restrictions, when I copy dir1 to pos0 with this angle?
        bool isRestricted0 = false;
        {
            isRestricted0 = allowedAngles.find(dir_gid1) != allowedAngles.end();
        }

        _Scalar angdiff0 = _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>( prim0, prim1,
                                                                                     /*isRestricted0 ? allowedAngles[dir_gid1] : */ angles,
                                                                                     &closest_angle_id0 );
        // prim0 needs to be close to prim1 in angle
        //if ( notPROMOTED(gid0,lid0) ) // don't bias small patches, copy best match to promoted patches
            add0 &= angdiff0 < angle_limit;

        // was this direction-angle pair already covered?
        {
            if ( copied[gid0].find(DidAid(dir_gid1,closest_angle_id0)) != copied[gid0].end() )
                add0 = false;
        }

        // are we still adding it?
        if ( !add0 )
            return false;

        // can candidate be created?
        _PrimitiveT cand0;
        if ( !prim0.generateFrom(cand0, prim1, closest_angle_id0, angles, _Scalar(1.)) )
            return false;

        // TEST triplets
        if ( tripletSafe ) // cancel candidate, if not ideally angled with any of the candidates with same dId already
        {
            for ( typename _PrimitiveContainerT::ConstIterator it(out_prims); it.hasNext(); it.step() )
            {
                if ( it.getDid() == dir_gid1 )
                {
                    //_Scalar tmpAng = GF2::angleInRad( cand0.template dir(), it->template dir() );
                    _Scalar tmpAng = _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>( cand0, *it,
                                                                                               angles );
                    if ( tmpAng < 0. ) { std::cout << "[" << __func__ << "]: " << "halt, tmpAng < 0. " << std::endl; throw std::runtime_error("asdf"); }
                    if ( tmpAng > _Scalar(1e-4) )
                    {
//                        std::cout << "triplet cancel: "
//                                  << "<" << gid0 << "," << dir_gid1 << "," << cand0.getTag(_PrimitiveT::TAGS::GEN_ANGLE) << "> "
//                                  << tmpAng * 180. / M_PI << " with "
//                                  << " <" << it.getGid() << "," << it.getDid() << "," << it->getTag(_PrimitiveT::TAGS::GEN_ANGLE) << ">"
//                                  << ", existing angle: " << _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>( prim1, *it,
//                                                                                                                       angles )
//                                  << "\n";
                        add0 = false;
                        break;
                    } //...if angle exceeds limits
                } //...if same dId
            } //...for all primitives

            // overwrite with a parallel copy of the best datafit one
            if ( !add0 )
            {
                typename _PrimitiveT::ExtentsT  extrema;
                std::vector<PidT>               population;
                processing::getPopulationOf( population, gid0, points );
                _PrimitiveT const*              bestPrim = NULL;
                Scalar                          bestDataCost = std::numeric_limits<Scalar>::max();

                for ( typename _PrimitiveContainerT::ConstIterator it(out_prims); it.hasNext(); it.step() )
                {
                    if ( it.getDid() == dir_gid1 )
                    {
                        // record this primitive with the same dId, and see, if it has best datacost.
                        // If the candidate get's cancelled, this might be the primitive,
                        // from which the direction we want to copy

                        _PrimitiveT tmp;
                        prim0.generateFrom( tmp, *it, 0, angles, _Scalar(1.) );
                        tmp.template getExtent<PointPrimitiveT>( extrema, points, scale, &population, true );
                        if ( !extrema.size() )
                            std::cout << "asdf" << std::endl;
                        else
                        {
                            Scalar distSum( 0. );
                            for ( int pid_id = 0; pid_id != population.size(); ++pid_id )
                                distSum += cand0.getFiniteDistance( extrema, points[ population[pid_id] ].template pos() );
                            if ( distSum < bestDataCost )
                            {
                                bestDataCost = distSum;
                                bestPrim     = &(*it);
                            }
                        } //...if extrema exists for dId
                    } //...if same dId
                } //...for all prims

                if ( bestPrim )
                {
                    closest_angle_id0 = 0; // parallel
                    if ( !prim0.generateFrom(cand0, *bestPrim, 0, angles, _Scalar(1.)) )
                        return false;
                    if ( copied[gid0].find(DidAid(dir_gid1,closest_angle_id0)) != copied[gid0].end() )
                        add0 = false;
                    else
                        add0 = true;
                } //...if bestPrim
            } //...if tripletCancel
        } //...tripletSafe

        // are we still adding it?
        if ( !add0 ) // changed by Aron on 9/04/15
            return false;

        // (1) if not restricted ? output, restrict
        // (2)                   : if angle allowed ? output
        // (3)                                      : if (generator not 0,180) and (not alias created) ? create alias

        // cache
        const _Scalar bestAngle = angles[closest_angle_id0];
        // guess generator
        AnglesT generators;
        genAngles( generators, bestAngle, angle_gens_in_rad );
        if ( generators.size() > 1 ) { throw new CandidateGeneratorException("[addCandidate] generators.size > 1, are you sure about this?"); }
        else if ( !generators.size() )
        {
            _Scalar genAngle = _PrimitiveT::GEN_ANGLE_VALUES::UNSET;
            if ( prim1.getTag(_PrimitiveT::TAGS::GEN_ANGLE) != _PrimitiveT::GEN_ANGLE_VALUES::UNSET )
                genAngle = prim1.getTag( _PrimitiveT::TAGS::GEN_ANGLE );
            else if ( allowedAngles.find(dir_gid1) != allowedAngles.end() )
            {
                //AnglesT _angleGens;
                //deduceGenerators( _angleGens, allowedAngles.at(dir_gid1) );
                genAngle = allowedAngles.at(dir_gid1).at(1);
            }

            if ( genAngle != _PrimitiveT::GEN_ANGLE_VALUES::UNSET )
            {
                cand0.setTag( _PrimitiveT::TAGS::GEN_ANGLE, genAngle );
                if ( verbose ) std::cout << "[" << __func__ << "]: " << "assuming gen_angle copy: " << cand0.getTag(_PrimitiveT::TAGS::GEN_ANGLE) << std::endl;
                generators.push_back( genAngle );
            } //...if genAngle filled
        } //...if generators.empty()

        bool added0   = false; // debug
        bool doOutput = false;

        // (1) If not restricted
        if ( !isRestricted0 )
        {
            // (1.1) Output
            doOutput = true;

            // (1.2) Restrict
            // if not 0 or 180, so there's a clear generator
            if ( generators.size() )
                // assign direction to generator
                angles::appendAnglesFromGenerators( /*        out: */ allowedAngles[dir_gid1]
                                                  , /* generators: */ generators
                                                  , /*   no_paral: */ false
                                                  , /*    verbose: */ false
                                                  , /*      inRad: */ true );
        }
        // (2) If restricted --> If angle found is in the restriction
        //else if ( std::find(allowedAngles[dir_gid1].begin(), allowedAngles[dir_gid1].end(), bestAngle) != allowedAngles[dir_gid1].end() )
        else if ( angles::findAngle( bestAngle, allowedAngles[dir_gid1]) )
        {
            // (2.1) Output
            doOutput = true;
        }
        // (3) Restricted, and not to the found angle --> If generator exists, and no alias exists yet, make one
        else if (    aliases && generators.size()                                                   // there exists a clear generator (i.e. 60 or 90, NOT 0, or 180) AND
                  && (    (  aliases          ->find(dir_gid1     ) == aliases            ->end())  // AND (there's no alias to this did yet
                         || (std::find_if( (*aliases)[dir_gid1].begin()                             //      OR the aliases don't have this generator yet)
                                         , (*aliases)[dir_gid1].end()
                                         , [&generators](typename _AliasesT::mapped_type::value_type const& elem){ return std::abs(elem.first-generators[0]) < halfDeg;} )
                             == (*aliases)[dir_gid1].end() )
                     )
                )
        {
            // (3.1)
            (*aliases)[dir_gid1][ generators[0] ] = aliases->createAlias( gid0, lid0, generators[0], prim0 );
        }

        // Output
        if ( doOutput )
        {
            if ( generators.size() > 1 ) { throw new CandidateGeneratorException("[addCandidate] generators.size > 1, are you sure about this?"); }
            added0 = output( cand0, out_prims, generated, copied, nLines, promoted, gid0, lid0, closest_angle_id0, generators.size() ? generators[0] : _PrimitiveT::GEN_ANGLE_VALUES::UNSET );
            if ( gid0 == 50 || gid1 == 50 )
                std::cout << "adding at " << gid0 << "," << dir_gid0 << " from " << gid1  << ", " << dir_gid1 << std::endl;
        }

        return added0;
    }

    template <class _PrimitiveT>
    struct AliasT
    {
        typedef typename _PrimitiveT::Scalar Scalar;

        AliasT() : _gid( _PrimitiveT::STATUS_VALUES::INVALID ), _lid( -1 ), _generatorAngle( -1 ), _prim( NULL ) {}

        AliasT( GidT const gid, GidT const lid, Scalar const generatorAngle, _PrimitiveT const& prim )
            : _gid( gid ), _lid( lid ), _generatorAngle( generatorAngle ), _prim( &prim ) {}
        AliasT( AliasT<_PrimitiveT> const& other )
            : _gid( other._gid ), _lid( other._lid ), _generatorAngle( other._generatorAngle ), _prim( other._prim ) {}

//        inline AliasT<_PrimitiveT> const& operator=( AliasT<_PrimitiveT> const& other )
//        {
//            if ( this == &other ) return *this;
//        }

        GidT    _gid, _lid;
        Scalar _generatorAngle; //!< \todo Is also stored in primitive now ( getTag( _PrimitiveT::ANGLE) )
        _PrimitiveT const* _prim;
    };

    template <class _PrimitiveT, typename _Scalar>
    struct AliasesT : public std::map< DidT, std::map<_Scalar,AliasT<_PrimitiveT> > >
    {
        inline AliasT<_PrimitiveT> createAlias( GidT const gid, GidT const lid, _Scalar generatorAngle, _PrimitiveT const& prim )
        {
            return AliasT<_PrimitiveT>( gid, lid, generatorAngle, prim );
        }
    };

    /*! \brief Stores candidate count for each direction ID. Is sortable by count.
     *         Used to filter out directions with only one instance amongst candidates.
     *         If there's only one candidate, we don't want it to be chosen,
     *         because it was not copied anywhere.
     * \note   Added by Aron on 7/1/2015 13:47
     */
    class CompareDirHistElements {
        public:
            typedef std::pair<DidT,LidT> ParentT;
            bool operator()(ParentT const& a, ParentT const& b ) { return a.second < b.second; }
    };

//    struct DirHistElemT : std::pair<DidT,LidT>
//    {
//        typedef std::pair<DidT,LidT> ParentT;
//        bool operator<( const ParentT& other ) const
//        {
//            return this->second < other.second;
//        }
//    };

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
    CandidateGenerator::generate( _PrimitiveContainerT                   & outPrims
                                , _PrimitiveContainerT            /*const*/& inPrims //non-const, becuase some of the primitives get promoted to large patches
                                , _PointContainerT                  const& points   // non-const to be able to add group tags
                                , _Scalar                           const  scale
                                , AnglesT                           const& angles
                                , CandidateGeneratorParams<_Scalar> const& params
                                , _Scalar                           const  smallThreshMult
                                , AnglesT                           const  angle_gens_in_rad
                                , bool                              const  safe_mode_arg
                                , int                               const  var_limit
                                , bool                              const  keepSingles
                                , bool                              const  allowPromoted
                                , bool                              const  tripletSafe
                                , bool const noAngleGuess
                                )
    {
        const bool verbose = true;

        // _________ typedefs _________

        typedef typename _PointContainerT::value_type                      _PointPrimitiveT;
        typedef typename _PrimitiveContainerT::const_iterator              outer_const_iterator;
        typedef typename _PrimitiveContainerT::iterator                    outer_iterator;
        typedef typename _PrimitiveContainerT::mapped_type                 InnerContainerT;
        //typedef typename outer_const_iterator::value_type::const_iterator inner_const_iterator;
        typedef typename InnerContainerT::const_iterator                   inner_const_iterator;
        typedef typename InnerContainerT::iterator                         inner_iterator;
        typedef          std::map< GidLid, LidT >                          GeneratedMapT;
        typedef          std::pair<std::pair<GidT,LidT>,unsigned long>     GeneratedEntryT;
        typedef          containers::PrimitiveContainer<_PrimitiveT>       PrimitiveMapT; //!< Iterable by single loop using PrimitiveMapT::Iterator

        // cache, so that it can be turned off
        bool safe_mode = safe_mode_arg;
        AliasesT<_PrimitiveT,_Scalar> aliases;

        // _________ error checks _________

        if ( outPrims.size() ) std::cerr << "[" << __func__ << "]: " << "warning, out_lines not empty!" << std::endl;
        if ( params.patch_population_limit < 0 ) { std::cerr << "[" << __func__ << "]: " << "error, popfilter is necessary!!!" << std::endl; return EXIT_FAILURE; }
        if ( angles[0] != _Scalar(0.) )
            throw new CandidateGeneratorException("angles[0] = 0 always! we need it for parallel");
        if ( params.small_mode != CandidateGeneratorParams<_Scalar>::SmallPatchesMode::IGNORE )
            throw new CandidateGeneratorException("small option is out of order, use IGNORE, meaning ignore small");

        // _________ variables _________

        LidT                                nlines = 0; // output size
        // filter already copied directions
        std::map< GidT, std::set<DidAid> >   copied; // [gid] = [ <dir0,angle_id0>, <dir0,angle_id1>, <dir2,angle_id0>, ... ]
        // modified angular similarity
        const _Scalar                       angle_limit( params.angle_limit / params.angle_limit_div );
        GeneratedMapT                       generated; // records, how many extra candidates the input primitive generated in the output


        // count patch populations
        std::cout << "[" << __func__ << "]: " << "populations start" << std::endl; fflush(stdout);
        GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        processing::getPopulations( populations, points );
        std::cout << "[" << __func__ << "]: " << "populations end" << std::endl; fflush(stdout);

        // _________ (1) promotion _________
        std::cout << "[" << __func__ << "]: " << "promotion start" << std::endl; fflush(stdout);
        // upgrade primitives to large, and copy large primitives to output
        // added on 21/09/2014 by Aron
        std::set<GidLid> promoted; // Contains primitives, that were small, but now are large (active)
        const _Scalar smallThresh = smallThreshMult * scale;
        {
            Eigen::Matrix<_Scalar,Eigen::Dynamic,1> spatialSignif(1,1); // cache variable

            // either all (patches), or none (second iteration) have to be unset
            int unset_count = 0, input_count = 0;
            // for all input primitives
            // #pragma omp parallel for private(spatialSignif)
            for ( outer_iterator outer_it0 = inPrims.begin(); outer_it0 != inPrims.end(); ++outer_it0 )
            {
                int gid = (*outer_it0).first;
                int lid = 0;
                for ( inner_iterator inner_it0 = (*outer_it0).second.begin(); inner_it0 != (*outer_it0).second.end(); ++inner_it0, ++lid )
                {
                    // cache vars
                    _PrimitiveT& prim         = *inner_it0;
                    int          prim_status  = prim.getTag(_PrimitiveT::TAGS::STATUS);

                    // bookkeeping
                    if ( prim_status == _PrimitiveT::STATUS_VALUES::UNSET )
                    {
//#                       pragma omp critical CRIT_UNSET_COUNT
                        ++unset_count;
                    }

                    // if unset (first iteration) or small (second+ iteration ), threshold it
                    if (    (prim_status == _PrimitiveT::STATUS_VALUES::SMALL)      // promotable OR
                         || (prim_status == _PrimitiveT::STATUS_VALUES::UNSET) )    // first ranking has to be done, since this is first iteration patchify output
                    {
                        // Promote: check, if patch is large by now
                        if (    ((populations.find(gid) != populations.end()) && (populations[gid].size() > params.patch_population_limit) )
                             && (prim.getSpatialSignificance( spatialSignif, points, scale, &(populations[gid]))(0) >= smallThresh  )   )
                        {
                            // store primitives, that have just been promoted to large from small
                            if ( (prim_status == _PrimitiveT::STATUS_VALUES::SMALL) )
                            {
//#                               pragma omp critical CRIT_UNSET_PROMOTED
                                promoted.insert( GidLid(gid,lid) );
                            }

                            prim.setTag( _PrimitiveT::TAGS::STATUS, _PrimitiveT::STATUS_VALUES::UNSET ); // set to unset, so that formulate can distinguish between promoted and
                        }
                        // Demote, if first iteration
                        else if (prim_status == _PrimitiveT::STATUS_VALUES::UNSET) // this should only happen in the first iteration
                        {
                            //std::cout << "demoting" << spatialSignif(0) << std::endl;
                            prim.setTag( _PrimitiveT::TAGS::STATUS, _PrimitiveT::STATUS_VALUES::SMALL );
                        }
                        else if (prim_status != _PrimitiveT::STATUS_VALUES::SMALL) // this should only happen in the first iteration
                        {
                            std::cout << "large patch( status: " << prim_status << "): " << spatialSignif(0) << std::endl;
                        }
//                        else
//                            std::cout << "not promoting size " << spatialSignif(0) << " < " << smallThresh << std::endl;

                        // update cache
                        prim_status = prim.getTag( _PrimitiveT::TAGS::STATUS );
                    } //...if not large

                    // copy to output
                    {
                        // copy input - to keep CHOSEN tag entries, so that we can start second iteration from a selection
                        auto &added = containers::add( outPrims, gid, *inner_it0 );

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
                            copied[ gid ].insert( DidAid(inner_it0->getTag(_PrimitiveT::TAGS::DIR_GID),0) );
                            //std::cout << "added to copied[" << gid << "]:" << inner_it0->getTag(_PrimitiveT::TAGS::DIR_GID) << std::endl;
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
            //if ( !promoted.size() && safe_mode )
            //{
            //    std::cout << "[" << __func__ << "]: " << "\n\nDISABLING safe mode, since no promoted\n" << std::endl;
            //    safe_mode = false;
            //}
        } //...promote
        std::cout << "[" << __func__ << "]: " << "populations end" << std::endl; fflush(stdout);

        // _________ (2) statistics _________
        std::cout << "[" << __func__ << "]: " << "angle stats start" << std::endl; fflush(stdout);
        // ANGLES
        std::map<DidT,AnglesT> allowedAngles;      // final storage to store allowed angles
        //AnglesT                angle_gens_in_rad;
        {
            std::map< _Scalar, AnglesT> angleSet;
            // ____ (1) store generators from earlier iterations in allowedAngles, so that they are not rediscovered ____
            for ( typename PrimitiveMapT::ConstIterator primIt(inPrims); primIt.hasNext(); primIt.step() )
            {
                if ( primIt->getTag( _PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::SMALL ) continue;
                // _PrimitiveT prim = *primIt;
                _Scalar const angleGen = primIt->getTag( _PrimitiveT::TAGS::GEN_ANGLE );

                // skip, if not meaningful
                if ( std::abs( angleGen - _PrimitiveT::GEN_ANGLE_VALUES::UNSET) < _Scalar(1.e-6) )
                    continue;

                DidT    const dId      = primIt->getTag( _PrimitiveT::TAGS::DIR_GID   );

                // only add, if has not been added before (gen_angle->did is stored at every primitive with that did redundantly for now)
                if ( allowedAngles.find(dId) == allowedAngles.end() )
                {
                    auto const angleSetIt = angleSet.find( angleGen );
                    if ( angleSetIt != angleSet.end() )
                    {
                        allowedAngles[dId] = angleSetIt->second;
                        continue;
                    }

                    // create a vector with this single element
                    AnglesT generators( {angleGen} );
                    // store angles from single generator for did
                    angles::appendAnglesFromGenerators( /*      out: */ allowedAngles[dId]
                                                      , /*       in: */ generators
                                                      , /* no_paral: */ false               // allow parallel
                                                      , /*  verbose: */ false
                                                      , /*    inRad: */ true );             // stored in rad!
                    angleSet[ angleGen ] = allowedAngles[dId];
                } //...if dId did not exist
            } //...add allowedAngles from primitive

            // ____ (2) Deduce rest from existing relations in the currently active landscape ____
            if ( !noAngleGuess )
            {
                // 2.1 Understand from angles, which generators we got from cli (i.e. 90, or 60,90)
                //deduceGenerators<_Scalar>( angle_gens_in_rad, angles );
                // log
                std::cout<<"[" << __func__ << "]: " << "angle_gens_in_rad:";for(size_t vi=0;vi!=angle_gens_in_rad.size();++vi)std::cout<<angle_gens_in_rad[vi]*180./M_PI<<" ";std::cout << "\n";

                // 2.2 estimate direction cluster angles
                {
                    DirAngleMapT dirAngles; // intermediate storage to collect dir angle distributions, <did, <angle_id, count> >
                    // create functor
                    DirAnglesFunctorOuter<_PrimitiveT,_PrimitiveContainerT,AnglesT> outerFunctor(
                                inPrims, angles, /* output reference */ dirAngles );

                    // run functor -> outputs { <did, [ <angle_gen,count>,... ]>, ... } to dirAngles
                    processing::filterPrimitives<_PrimitiveT, InnerContainerT/*typename _PrimitiveContainerT::mapped_type*/>(
                                inPrims, outerFunctor );

                    // dirAngles -> allowedAngles { <did, [ 0, angle_i, ..., M_PI ] >, ... }
                    selectAngles( /* out: */ allowedAngles, /* in: */ dirAngles, angles, angle_gens_in_rad );
                }
            }
        }
        std::cout << "[" << __func__ << "]: " << "angle stats stop" << std::endl; fflush(stdout);

        /* Status of primitives at this stage are:
         * ACTIVE      for large primitives to use
         * TAG_UNSET   for promoted primitives to give directions to
         * SMALL       for primitives to ignore
         */

        // _________ (3) generation _________
        std::cout << "[" << __func__ << "]: " << "generate start" << std::endl; fflush(stdout);

        DidT maxDid = 0; // collects currently existing maximum cluster id (!small, active, all!)

        GidT gid0, gid1, dir_gid0, dir_gid1, lid0, lid1;
        gid0 = gid1 = dir_gid0 = dir_gid1 = _PrimitiveT::TAG_UNSET; // group tags cached
        for ( outer_const_iterator outer_it0  = inPrims.begin(); outer_it0 != inPrims.end(); ++outer_it0 )
        {
            gid0 = -2; // -1 is unset, -2 is unread
            lid0 =  0;
            for ( inner_const_iterator inner_it0  = (*outer_it0).second.begin(); inner_it0 != (*outer_it0).second.end(); ++inner_it0, ++lid0 )
            {
                // cache outer primitive
                _PrimitiveT const& prim0 = *inner_it0;
                dir_gid0                 = prim0.getTag( _PrimitiveT::TAGS::DIR_GID );

                if ( dir_gid0 > maxDid ) maxDid = dir_gid0; // small, active, uninited!

                // cache group id of patch at first member
                     if ( gid0 == -2               )  gid0 = prim0.getTag( _PrimitiveT::TAGS::GID ); // store gid of first member in patch
                else if ( gid0 != outer_it0->first )  std::cerr << "[" << __func__ << "]: " << "Not good, prims under one gid don't have same GID..." << std::endl;

                for ( outer_const_iterator outer_it1  = outer_it0; outer_it1 != inPrims.end(); ++outer_it1 )
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

                        // cache group id of patch at first member
                             if ( gid1 == -2               )    gid1 = prim1.getTag( _PrimitiveT::TAGS::GID );
                        else if ( gid1 != outer_it1->first )    std::cerr << "[" << __func__ << "]: " << "Not good, prims under one gid don't have same GID in inner loop..." << std::endl;

                        addCandidate<_PrimitivePrimitiveAngleFunctorT>(
                                    prim0, prim1, lid0, lid1, safe_mode, allowPromoted, angle_limit, angles, angle_gens_in_rad, promoted,
                                    allowedAngles, copied, generated, nlines, outPrims, points, scale, &aliases, tripletSafe );
                        addCandidate<_PrimitivePrimitiveAngleFunctorT>(
                                    prim1, prim0, lid1, lid0, safe_mode, allowPromoted, angle_limit, angles, angle_gens_in_rad, promoted,
                                    allowedAngles, copied, generated, nlines, outPrims, points, scale, &aliases, tripletSafe );
                    } //...for l3
                } //...for l2
            } //...for l1
        } //...for l0
        std::cout << "[" << __func__ << "]: " << "generate end" << std::endl; fflush(stdout);

        // ___________ (4) ALIASES _______________
        std::cout << "[" << __func__ << "]: " << "alias start" << std::endl; fflush(stdout);

        // [did][angle] = AliasT( gid, lid, prim )
        for ( auto aliasIt = aliases.begin(); aliasIt != aliases.end(); ++aliasIt )
        {
            // first: dir_gid
            // second: std::map<angle, AliasT>
            for ( auto angleIt = aliasIt->second.begin(); angleIt != aliasIt->second.end(); ++angleIt )
            {
                // first: angle
                // second: AliasT

                // copy original primitive
                _PrimitiveT prim0 = *( angleIt->second._prim );
                const DidT did = ++maxDid;
                prim0.setTag( _PrimitiveT::TAGS::DIR_GID, did );
                const GidT gid0 = angleIt->second._gid;
                const LidT lid0 = angleIt->second._lid;

                if ( !output( prim0, outPrims, generated, copied, nlines, promoted, gid0, lid0, 0, prim0.getTag(_PrimitiveT::TAGS::GEN_ANGLE)) )
                {
                    std::cerr << "Alias could not be added" << std::endl;
                    throw new CandidateGeneratorException("Alias could not be added");
                } //...if could not output
                else
                    if ( verbose ) std::cout << "\nAdded alias to " << angleIt->second._prim->getTag(_PrimitiveT::TAGS::DIR_GID)
                              << " as " << outPrims[gid0].back().toString()
                              << " <" << outPrims[gid0].back().getTag( _PrimitiveT::TAGS::GID ) << ","
                              << outPrims[gid0].back().getTag( _PrimitiveT::TAGS::DIR_GID ) << ">"
                              << std::endl;

                // record allowed alias
                AnglesT tmpGenerators( {angleIt->first} );
                angles::appendAnglesFromGenerators( /*        out: */ allowedAngles[ did ]
                                                  , /* generators: */ tmpGenerators
                                                  , /*   no_paral: */ false
                                                  , /*    verbose: */ false
                                                  , /*      inRad: */ true );

                // add all allowed copies
                for ( outer_const_iterator outer_it  = inPrims.begin(); outer_it != inPrims.end(); ++outer_it )
                {
                    lid1 =  0;
                    for ( inner_const_iterator inner_it = (*outer_it).second.begin(); inner_it != (*outer_it).second.end(); ++inner_it, ++lid1 )
                    {
                        // cache outer primitive
                        _PrimitiveT const& prim1 = *inner_it;

                        // copy prim0 (the alias) to all compatible receivers given allowedAngles.
                        addCandidate<_PrimitivePrimitiveAngleFunctorT,AliasesT<_PrimitiveT,_Scalar> >(
                            prim1, prim0, lid1, lid0, safe_mode, allowPromoted, angle_limit, angles, angle_gens_in_rad, promoted,
                            allowedAngles, copied, generated, nlines, outPrims, points, scale, nullptr, tripletSafe );

                    } //...inner for
                } //...outer for
            } //...for all angles
        } //...for all aliases
        std::cout << "[" << __func__ << "]: " << "alias end" << std::endl; fflush(stdout);

        // ___________ (5) LIMIT VARIABLES _______________
        std::cout << "[" << __func__ << "]: " << "limit vars start" << std::endl; fflush(stdout);
        int ret = EXIT_SUCCESS;
        if ( (var_limit > 0) && (nlines > var_limit) )
        {
            LidT active_count = nlines; // records how many variables actives generated among themselves, these we cannot filter

            // (1) sort promoted patches descending based on how many variables they generated
            // in: generated, containing < <gid,lid> , generated count> entries for promoted locations only (no actives)
            // out: ranks, containing a sorted list of generated
            std::vector< GeneratedEntryT > ranks;
            for ( auto it = generated.begin(); it != generated.end(); ++it )
            {
                // in_lines[it.first.first][it.first.second].getTag( _PrimitiveT::TAGS::STATUS ) )
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
                ret = -1;
                // demote!
                LidT demoted = 0;
                for ( auto pit = promoted.begin(); pit != promoted.end(); ++pit )
                {
                    _PrimitiveT &prim = inPrims.at( (*pit).first  )
                                               .at( (*pit).second );

                    if ( prim.getTag(_PrimitiveT::TAGS::STATUS ) == _PrimitiveT::STATUS_VALUES::ACTIVE )
                        std::cerr << "demoting active " << prim.getTag(_PrimitiveT::TAGS::GID) << "(gid) ????" << std::endl;
                    else
                        ++demoted;

                    prim.setTag(_PrimitiveT::TAGS::STATUS, _PrimitiveT::STATUS_VALUES::SMALL);
                }
                std::cout << "[" << __func__ << "]: " << "demoted " << demoted << " primitives" << std::endl;
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
                    LidT vars = active_count;
                    for ( LidT i = 0; i != ranks.size(); ++i )
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
                auto tmp = outPrims; outPrims.clear();
                cgen::FilterChosenFunctor<_PrimitiveT, _PrimitiveContainerT> functor( outPrims, chosen, promoted );
                nlines = processing::filterPrimitives<_PrimitiveT, InnerContainerT>( tmp, functor );
                // return remaining variable count
                ret = ranks.size() - chosen.size();
            } // if too many actives already
        } // filter
        std::cout << "[" << __func__ << "]: " << "limit vars end" << std::endl; fflush(stdout);

        // ___________ (6) Make sure allowed angles stick _______________
        std::cout << "[" << __func__ << "]: " << "allowed start" << std::endl; fflush(stdout);
        std::map<DidT,LidT> directionPopulation; // counts, how many candidates have a direction
        for ( typename PrimitiveMapT::Iterator primIt(outPrims); primIt.hasNext() && (ret != -1); primIt.step() )
        {
            // _PrimitiveT prim = *primIt;
            DidT    const dId      = primIt->getTag( _PrimitiveT::TAGS::DIR_GID   );
            _Scalar const angleGen = primIt->getTag( _PrimitiveT::TAGS::GEN_ANGLE );

            // only add, if has not been added before (gen_angle->did is stored at every primitive with that did redundantly for now)
            if ( allowedAngles.find(dId) != allowedAngles.end() && (angleGen == _PrimitiveT::STATUS_VALUES::UNSET) )
            {
//                std::cout << "[" << __func__ << "]: " << "primitive(" << primIt.getGid() << "," << dId << ") "
//                             << "does not have GEN_ANGLE stored: " << angleGen
//                             << ", but allowedAngles has entry [1]: " << allowedAngles[dId][1] << std::endl;
                primIt->setTag( _PrimitiveT::TAGS::GEN_ANGLE, allowedAngles[dId][1] );
            }

            // 6.2 Histogram direction ID-s, and list the ones, that are alone (will not get chosen anyway)
            directionPopulation[ dId ]++;
        } //...for outPrims

        // 6.3 Filter primitives with only one direction among candidates
        if ( !keepSingles  && (ret != -1) )
        {
            std::map<GidT,LidT> thrownCnt;
            auto tmp = outPrims; outPrims.clear();
            for ( typename PrimitiveMapT::Iterator primIt(tmp); primIt.hasNext(); primIt.step() )
            {
                // _PrimitiveT prim = *primIt;
                DidT const dId      = primIt->getTag( _PrimitiveT::TAGS::DIR_GID );
                if (    (directionPopulation[dId] > 1)
                     || (primIt->getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::SMALL)
                     || (    (populations.find(primIt.getGid()) != populations.end()      )
                          && (populations[primIt.getGid()].size() > params.patch_population_limit * 2) ) // add if, at least double poplimit points in it
                   )
                {
                    containers::add( outPrims, primIt.getGid(), *primIt );
                }
                else
                {
                    thrownCnt[primIt.getGid()]++;
                }
            }

            for ( auto it = thrownCnt.begin(); it != thrownCnt.end(); ++it )
            {
                if ( outPrims[ it->first ].size() == 0 )
                {
                    std::cerr << "[" << __func__ << "]: " << "throwing away everything from gid " << it->first << std::endl;
                    //throw std::runtime_error("gid");
                }
            }
        }
        std::cout << "[" << __func__ << "]: " << "allowed end" << std::endl; fflush(stdout);

        // log
        std::cout << "[" << __func__ << "]: " << "finished generating, we now have " << nlines << " candidates" << std::endl;

        return ret;
    } // ...CandidateGenerator::generate()

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
    int
    CandidateGenerator::generateCli( int    argc
                                   , char** argv )
    {
        typedef typename _PrimitiveContainerT::value_type InnerPrimitiveContainerT;
        //typedef typename PointContainerT::value_type PointPrimitiveT;
        int err = EXIT_SUCCESS;

        CandidateGeneratorParams<_Scalar> generatorParams;
        std::string                 cloud_path              = "./cloud.ply";
        AnglesT                     angleGens( {_Scalar(90.)} );
        std::string                 mode_string             = "representative_sqr";
        std::vector<std::string>    mode_opts               = { "representative_sqr" };
        std::string                 input_prims_path        = "patches.csv";
        std::string                 associations_path       = "points_primitives.csv";
        bool                        keepSingles             = false; // throw away single dId-s
        bool                        allowPromoted           = false;
        bool                        tripletSafe             = false;

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
            pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angleGens );
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
                generatorParams.small_mode = static_cast<typename CandidateGeneratorParams<_Scalar>::SmallPatchesMode>( small_mode );
            }

            // keepSingles
            {
                keepSingles = pcl::console::find_switch( argc, argv, "--keep-singles" );
                allowPromoted= pcl::console::find_switch( argc, argv, "--allow-promoted" );
                tripletSafe = pcl::console::find_switch( argc, argv, "--triplet-safe");
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
                    std::cout << "\t [--angle-gens "; for(size_t vi=0;vi!=angleGens.size();++vi)std::cout<<angleGens[vi]<<","; std::cout << "]\n";
                    std::cout << "\t [--patch-pop-limit " << generatorParams.patch_population_limit << "]\n";
                    std::cout << "\t [--small-mode " << generatorParams.small_mode << "\t | 0: IGNORE, 1: RECEIVE_SIMILAR, 2: RECEIVE_ALL]\n";
                    std::cout << "\t [--no-paral]\n";
                    std::cout << "\t [--safe-mode]\n";
                    std::cout << "\t [--var-limit " << generatorParams.var_limit << "\t Decides how many variables we want as output. 0 means unlimited.]\n";
                    std::cout << "\t [--keep-singles " << (keepSingles?"YES":"NO") << "\t Decides, if we should throw away single directions]\n";
                    std::cout << "\t [--allow-promoted " << (allowPromoted?"YES":"NO") << "\t Decides, if we should allow promoted patches to distribute their directions]\n";
                    std::cout << "\t [--triplet-safe " << (tripletSafe?"YES":"NO") << "]\t Ensure, that perfect angle is respected for every member of DiD\n";
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
        AnglesT angleGensInRad;
        if ( EXIT_SUCCESS == err )
        {
            angles::appendAnglesFromGenerators( generatorParams.angles, angleGens, no_paral, false );
            for ( typename AnglesT::const_iterator it = angleGens.begin(); it != angleGens.end(); ++it )
                angleGensInRad.push_back( *it * M_PI / _Scalar(180.) );
        } //...read angles

        // Read points
        _PointContainerT points;
        if ( EXIT_SUCCESS == err )
        {
            err = io::readPoints<_PointPrimitiveT>( points, cloud_path );
            if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
        } //...read points

        std::vector<std::pair<GidT,LidT> > points_primitives;
        io::readAssociations( points_primitives, associations_path, NULL );
        for ( size_t i = 0; i != points.size(); ++i )
        {
            // store association in point
            points[i].template setTag( _PointPrimitiveT::TAGS::GID, points_primitives[i].first );
        }

        // read primitives
        //typedef std::map<GidT, InnerPrimitiveContainerT> PrimitiveMapT; // redefined by Aron 14.23 13/1/2015
        typedef containers::PrimitiveContainer<_PrimitiveT> PrimitiveMapT;
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
        int attempts = 0, attemptLimit = 1;
        if ( EXIT_SUCCESS == err )
        {
            // runs once, 0<1, unless too many actives in output,
            // in which case, runs twice (0<1, 1<2) with safe_mode on second run
            while ( attempts < attemptLimit )
            {
                ++attempts;

                ret = CandidateGenerator::generate< MyPrimitivePrimitiveAngleFunctor, MyPointPrimitiveDistanceFunctor, _PrimitiveT >
                        (   /* out: */ primitives
                          , /*  in: */ patches
                          , points
                          , generatorParams.scale
                          , generatorParams.angles
                          , generatorParams
                          , generatorParams.small_thresh_mult
                          , angleGensInRad
                          , generatorParams.safe_mode
                          , generatorParams.var_limit
                          , keepSingles
                          , allowPromoted
                          , tripletSafe
                          );

                //if ( err != EXIT_SUCCESS ) std::cerr << "[" << __func__ << "]: " << "generate exited with error! Code: " << err << std::endl;
                if ( ret > 0 )
                    std::cout << "not all patches were promoted( " << ret << " left), will need rerun on same threshold..." << std::endl;
                else if ( ret < 0 )
                {
                    std::cout << "rerunning in safe mode, all active..." << std::endl;
                    attemptLimit = 2;
                    generatorParams.safe_mode = 1;
                    primitives.clear();

                    // demote back all promoted
                    for ( typename containers::PrimitiveContainer<_PrimitiveT>::Iterator it(patches); it.hasNext(); it.step() )
                        if ( it->getTag(_PrimitiveT::TAGS::STATUS) == _PrimitiveT::STATUS_VALUES::UNSET )
                             it->setTag( _PrimitiveT::TAGS::STATUS, _PrimitiveT::STATUS_VALUES::SMALL );
                }
            }
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

//        if ( attempts > 1 )
//            return 1;
        return std::max(ret,err);
    } // ...CandidateGenerator::generateCli()

} // ... ns GF2

#undef isSMALL
#undef notSMALL
#undef isPROMOTED
#undef notPROMOTED

#endif // GF2_CANDIDATEGENERATOR_HPP
