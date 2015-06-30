#ifndef RAPTER_SUBSAMPLEPRIMITIVES_HPP
#define RAPTER_SUBSAMPLEPRIMITIVES_HPP

#include "rapter/simpleTypes.h"
#include "rapter/io/inputParser.hpp"
#include "rapter/processing/util.hpp"
#include "rapter/util/containers.hpp"
#include "rapter/util/impl/randUtil.hpp"
#include <chrono>

namespace rapter
{
    // Usage: .../Release/bin/toGlobFit --subsample-primitives 0.1 --pop-limit 100 --prims segments.csv --cloud cloud.ply -a points_segments.csv --scale 0.005
    template < class _PrimitiveVectorT
             , class _PrimitiveMapT
             , class _PointContainerT
             , class _PclCloudT>
    int subsamplePrimitives( int argc, char** argv )
    {
        using namespace rapter;

        typedef typename _PrimitiveMapT::mapped_type    InnerContainerT;
        typedef typename _PclCloudT::Ptr                PclPtrT;
        typedef typename InnerContainerT::value_type    PrimitiveT;
        typedef typename PrimitiveT::Scalar             Scalar;
        typedef typename _PointContainerT::PrimitiveT   PointPrimitiveT;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        srand( time(NULL) );

        _PointContainerT    points;
        PclPtrT             pclCloud;
        _PrimitiveVectorT   primitivesVector;
        _PrimitiveMapT      primitives;
        struct MyParams { Scalar scale; } params;

        rapter::parseInput<InnerContainerT,_PclCloudT>( points, pclCloud, primitivesVector, primitives, params, argc, argv, /* readAssoc: */ true );

        Scalar ratio( 1. );
        if ( rapter::console::parse_argument(argc,argv,"--subsample-primitives",ratio) < 0 || (ratio < FLT_EPSILON) )
        {
            std::cerr << "[" << __func__ << "]: " << "you need to provide subsample ratio after --subsample-primitives\n";
            std::cout << "example: " << "--subsample-primitives 0.1 --prim-limit 200 --prim-random-ratio 0.5 --pop-limit 100 --prims segments.csv --cloud cloud.ply -a points_segments.csv --scale 0.005" << std::endl;
            return EXIT_FAILURE;
        }
        int popLimit( 20 ); // don't cut planes less than this
        rapter::console::parse_argument( argc,argv,"--pop-limit", popLimit );
        std::cout << "[" << __func__ << "]: " << "keeping at least " << popLimit << " points for each plane, change with \"--pop-limit k\"" << std::endl;

        size_t primLimit( 0 ); // keep the first primLimit planes
        rapter::console::parse_argument( argc,argv,"--prim-limit", primLimit );
        Scalar primRandRatio( 0. ); // randomize the selection of primitives
        rapter::console::parse_argument( argc,argv,"--prim-random-ratio", primRandRatio );
        std::cout << "[" << __func__ << "]: " << "keeping " << primLimit << " planes, change with \"--prim-limit k\"" << std::endl;
        std::cout << "[" << __func__ << "]: " << primRandRatio * 100.0 << "% will be random (--prim-random-ratio f), the rest picked in descending population (point count) order\n";

        // fetch populations
        rapter::GidPidVectorMap populations;
        processing::getPopulations( populations, points );

        //typedef typename GF2::GidPidVectorMap::value_type       Gid_PidVector_Pair;
        typedef typename rapter::GidPidVectorMap::const_iterator   PopConstIterator;
        typedef          std::pair<PidT, GidT>                  MultiMapPair;

        /// sort by population
        std::multimap< PidT, GidT > gidsBySize;         // key: #points in GId, value: { GId => [ pid0, pid1, ...] }
        PidT                        sumCanBeCut ( 0 );  // sum of #points that belong to planes, that >popLimit
        for ( PopConstIterator it = populations.begin(); it != populations.end(); ++it )
        {
            if ( it->first == PrimitiveT::LONG_VALUES::UNSET )
                continue;
            // read population
            const PidT pop = it->second.size();

            // record
            gidsBySize.insert( MultiMapPair(pop, it->first) );

            // count, how many points we plan to cut from (not all, planes < popLimit we don't touch)
            if ( pop > popLimit )
                sumCanBeCut += pop;
        } //...for populations


        /// estimate per-plane-reduction
        std::map<GidT,PidT> cutHowMuch;                             // how much each Gid has to be cut
        std::set<GidT>      keepPrims;                              // which primitives are large enough to keep. Only used, if primLimit > 0
        const PidT          targetCut( std::floor(Scalar(1.-ratio) * sumCanBeCut) ); // how many points we were instructed to cut
        PidT                sumCut   ( 0 );                         // how many points we were plan       to cut (sanity check ~= targetCut)

        // pick (1. - primRandRatio) primitives from the top, and the rest randomly
        std::set<MultiMapPair> chosenGids;
        {
            ULidT largePrimDemand = std::floor( primLimit * (1. - primRandRatio) );
            std::cout << "largePrimDemand " << largePrimDemand << ", primlimiT: " << primLimit << std::endl;

            // iterate over all once
            auto it = gidsBySize.rbegin();
            // pick first k right away
            for ( ; it != gidsBySize.rend() && (chosenGids.size() < largePrimDemand); ++it )
                chosenGids.insert( *it );

            // enlist the rest to pick randomly from
            std::vector<MultiMapPair> gidMix; // contains the "small" primitives, to randomly sample from
            for ( ; it != gidsBySize.rend(); ++it )
                gidMix.push_back( *it );
            std::cout << "gidMix.size(): " << gidMix.size() << ", " << chosenGids.size() << std::endl;

            // shuffle the "small" primitives
            shuffle( gidMix.begin(), gidMix.end(), std::default_random_engine(seed) );

            // pick the first (primLimit-k) to fill chosen Gids
            for ( auto gidMixIt = gidMix.begin(); gidMixIt != gidMix.end() && (chosenGids.size() < primLimit); ++gidMixIt )
            {
                chosenGids.insert( *gidMixIt );
                std::cout << "chosenGids.size(): " << chosenGids.size() << std::endl;
            }

        }

        // actually pick them, and subsample their population
        LidT primId( 0 );
        //for ( auto it = gidsBySize.rbegin(); it != gidsBySize.rend(); ++it, ++primId )
        for ( auto it = chosenGids.begin(); it != chosenGids.end(); ++it, ++primId )
        {
            const PidT population = it->first;
            const GidT gid        = it->second;
            //std::cout << "[#points:" << pop << "]: <gid," << gid << ">, #points:" << populations.at( it->second ).size() << std::endl;

            // remove primitive, if
            //  a. there's a primitive limit AND
            //  b. already exceeded that limit OR
            //  c. this is large primitive we randomly don't keep and leave space for a random one later: rand < primRandRatio
            if ( primLimit && (keepPrims.size() >= primLimit) )
            {
                std::cerr << "[" << __func__ << "]: " << " too many primitives chosen..." << std::endl;
                continue;
            }
            else
            {
                keepPrims.insert( gid );
                std::cout << "keeping " << gid << "/" << gidsBySize.size() << ", keepPrims.size(): " << keepPrims.size() << std::endl;
            }

            // don't cut points from too small planes
            if ( population < popLimit )
                continue;

            // this planes contribution to the cut
            Scalar  localRatio( Scalar(population) / Scalar(sumCanBeCut) );
            PidT    localN    ( std::floor(localRatio * targetCut) );
            if ( population - localN < popLimit )
                localN = population - popLimit;

            if ( localN && (population - localN < popLimit) )
                std::cerr << "[" << __func__ << "]: " << "warning, cutting a plane to become too small: "
                          << population
                          << " --> " << population - localN
                          << " < " << popLimit
                          << std::endl;
            // record
            if ( localN )
            {
                if ( cutHowMuch.find(gid) != cutHowMuch.end() )
                    std::cerr << "[" << __func__ << "]: " << "warning, probably more planes in same Gid" << std::endl;
                cutHowMuch[ gid ] = localN;
            }

            // sanity check
            sumCut += localN;

            // log
            //std::cout << "localRat = " << population << "/" << sumCanBeCut << ": " << localRatio << "\t=> * " << targetCut << " = " << localN << std::endl;
        } //...for planes to cut

        // log
        std::cout << "[" << __func__ << "]: " << "cut " << sumCut << "/" << points.size() << " assignments == " << Scalar(sumCut)/Scalar(sumCanBeCut) * 100. << "%" << std::endl;

        /// do cut
        _PrimitiveMapT out;
        for ( typename _PrimitiveMapT::Iterator it(primitives); it.hasNext(); it.step() )
        {
            const GidT gid = it.getGid();

            auto popIt = populations.find( gid );
            if ( popIt == populations.end() )
            {
                std::cerr << "[" << __func__ << "]: " << " this shouldn't happen, population not found for gid " << gid << std::endl;
                continue;
            }
            auto population = popIt->second; // copy, because we reshuffle it later

            // (1) clear all point assignments, we delete the primitive
            if ( primLimit && (keepPrims.find(gid) == keepPrims.end()) )
            {
                for ( auto pidIt = population.rbegin(); pidIt != population.rend() ; ++pidIt )
                {
                    if ( points.at( /* pid: */ *pidIt ).getTag( PointPrimitiveT::TAGS::GID ) != gid )
                        std::cerr << "[" << __func__ << "]: " << "Point " << *pidIt << " is not actually assigned to plane (gid " << gid << "), should not happen" << std::endl;
                    points.at( /* pid: */ *pidIt ).setTag( PointPrimitiveT::TAGS::GID, PointPrimitiveT::LONG_VALUES::UNSET );
                }
                continue;
            }

            // store primitive
            containers::add( out, gid, *it );
            std::cout << "adding " << gid << ", now " << out.size() << std::endl;

            // (2) just copy primitive to output, no assigned points to remove
            auto const cutIt = cutHowMuch.find( gid );
            if ( cutIt == cutHowMuch.end() )
            {
                continue;
            }

            // (3) keep, but remove some of the assigned points
            const PidT cut = cutIt->second;

            // estimate percentage of points that need to be cut from all of the points assigned to this primitive
            Scalar  cutRatio( Scalar(cut) / Scalar( population.size()) );
            // number of points to cut
            PidT    toCut   ( std::ceil(cutRatio * Scalar(population.size())) );

            // log/debug
            if ( cutRatio < ratio * 0.01  || cutRatio > 1.01 ) std::cerr << "cutratio: " << cutRatio << std::endl;

            // reorder assigned points to randomize
            shuffle ( population.begin(), population.end(), std::default_random_engine(seed) );
            // keep track of cut points
            PidT cnt   = 0;
            // iterate through shuffled, assigned points
            for ( auto pidIt = population.rbegin(); pidIt != population.rend() && cnt != toCut; ++pidIt, ++cnt )
            {
                if ( points.at( /* pid: */ *pidIt ).getTag( PointPrimitiveT::TAGS::GID ) != gid )
                    std::cerr << "[" << __func__ << "]: " << "Point " << *pidIt << " is not actually assigned to plane (gid " << gid << "), should not happen" << std::endl;
                points.at( /* pid: */ *pidIt ).setTag( PointPrimitiveT::TAGS::GID, PointPrimitiveT::LONG_VALUES::UNSET );
            }
        } //...iterate primitives

        // check
        rapter::GidPidVectorMap populations2;
        processing::getPopulations( populations2, points );
        PidT                 sum(0);
        for ( auto it = populations2.begin(); it != populations2.end(); ++it )
        {
            if ( it->first == PrimitiveT::LONG_VALUES::UNSET )
                continue;

            sum+= it->second.size();
        }
        std::cout << "[" << __func__ << "]: "
                  << "Reduction: from " << points.size() << " points, and " << sumCanBeCut << " assignments --> "
                  << sum << " remaining assignments  == "
                  << Scalar(sum) / points.size() * 100. << "% points are active in the output, "
                  << Scalar(sum) / sumCanBeCut   * 100. << "% assignments remained" << std::endl;

        for ( auto it = gidsBySize.begin(); it != gidsBySize.end(); ++it )
        {
            const PidT population = it->first;
            //const GidT gid = it->second;
            if ( population < popLimit )
                continue;
            //std::cout << "oldPop[" << gid << "]: " << populations[gid].size() << ", new: " << populations2[gid].size() << std::endl;
        }

        // save
        std::string input_prims_path = parsePrimitivesPath(argc, argv);
        std::string outName = boost::filesystem::path( input_prims_path).stem().string();
        std::stringstream ssPrims; ssPrims << outName << ".sub_" << ratio << "_" << primLimit << ".csv";
        rapter::io::savePrimitives<PrimitiveT,typename InnerContainerT::const_iterator>( out, ssPrims.str() );

        std::stringstream ssAssoc; ssAssoc << "points_" << outName << ".sub_" << ratio << "_" << primLimit << ".csv";
        rapter::io::writeAssociations<PointPrimitiveT>( points, ssAssoc.str() );

        std::cout << "[" << __func__ << "]: " << "results written to " << ssPrims.str() << "(#" << out.size() << ") and " << ssAssoc.str() << "\n";
        std::cout << "../runGlobfit.py -s " << params.scale << " -p " << ssPrims.str() << " -a " << ssAssoc.str() << "\n";

        return EXIT_SUCCESS;
    } //...subsamplePrimitives()
} //...ns rapter

#endif // RAPTER_SUBSAMPLEPRIMITIVES_HPP
