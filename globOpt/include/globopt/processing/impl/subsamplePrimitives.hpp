#ifndef GO_SUBSAMPLEPRIMITIVES_HPP
#define GO_SUBSAMPLEPRIMITIVES_HPP

#include "globfit2/simple_types.h"
#include "globfit2/io/inputParser.hpp"
#include "globfit2/processing/util.hpp"
#include "globfit2/util/containers.hpp"
#include <chrono>

namespace globopt
{
    template <typename _Scalar>
    inline _Scalar randf() { return rand() / _Scalar(RAND_MAX); }

    // Usage: .../Release/bin/toGlobFit --subsample-primitives 0.1 --pop-limit 100 --prims segments.csv --cloud cloud.ply -a points_segments.csv --scale 0.005
    template < class _PrimitiveVectorT
             , class _PrimitiveMapT
             , class _PointContainerT
             , class _PclCloudT>
    int subsamplePrimitives( int argc, char** argv )
    {
        using namespace GF2;

        typedef typename _PrimitiveMapT::mapped_type InnerContainerT;
        typedef typename _PclCloudT::Ptr             PclPtrT;
        typedef typename InnerContainerT::value_type PrimitiveT;
        typedef typename PrimitiveT::Scalar Scalar;
        typedef typename _PointContainerT::PrimitiveT PointPrimitiveT;

        srand( time(NULL) );
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        _PointContainerT    points;
        PclPtrT             pclCloud;
        _PrimitiveVectorT   primitivesVector;
        _PrimitiveMapT      primitives;
        struct MyParams { Scalar scale; } params;

        GF2::parseInput<InnerContainerT,_PclCloudT>( points, pclCloud, primitivesVector, primitives, params, argc, argv, /* readAssoc: */ true );

        Scalar ratio( 1. );
        if ( GF2::console::parse_argument(argc,argv,"--subsample-primitives",ratio) < 0 || (ratio < FLT_EPSILON) )
        {
            std::cerr << "[" << __func__ << "]: " << "you need to provide subsample ratio after --subsample-primitives\n";
            std::cout << "example: " << "--subsample-primitives 0.1 --pop-limit 100 --prims segments.csv --cloud cloud.ply -a points_segments.csv --scale 0.005" << std::endl;
            return EXIT_FAILURE;
        }
        int popLimit( 20 ); // don't cut planes less than this
        GF2::console::parse_argument( argc,argv,"--pop-limit", popLimit );
        std::cout << "keeping at least " << popLimit << " points for each plane, change with \"--pop-limit k\"" << std::endl;

        int primLimit( 0 ); // keep the first primLimit planes
        GF2::console::parse_argument( argc,argv,"--prim-limit", primLimit );
        std::cout << "keeping first " << primLimit << " planes, change with \"--prim-limit k\"" << std::endl;

        // fetch populations
        GF2::GidPidVectorMap populations;
        processing::getPopulations( populations, points );

        //typedef typename GF2::GidPidVectorMap::value_type       Gid_PidVector_Pair;
        typedef typename GF2::GidPidVectorMap::const_iterator   PopConstIterator;
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

        std::cout << "actualLim: " << popLimit << std::endl;
        for ( auto it = gidsBySize.rbegin(); it != gidsBySize.rend(); ++it )
        {
            const PidT population = it->first;
            const GidT gid = it->second;
            //std::cout << "[#points:" << pop << "]: <gid," << gid << ">, #points:" << populations.at( it->second ).size() << std::endl;

            if ( primLimit && (keepPrims.size() >= primLimit) )
                continue;
            else
                keepPrims.insert( gid );

            if ( population < popLimit )
                continue;

            // this planes contribution to the cut
            Scalar  localRatio( Scalar(population) / Scalar(sumCanBeCut) );
            PidT    localN    ( std::floor(localRatio * targetCut) );
            if (population - localN < popLimit)
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
            std::cout << "localRat = " << population << "/" << sumCanBeCut << ": " << localRatio << "\t=> * " << targetCut << " = " << localN << std::endl;
        } //...for planes to cut

        // log
        std::cout << "cut " << sumCut << "/" << points.size() << " == " << Scalar(sumCut)/Scalar(points.size()) << std::endl;

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
        GF2::GidPidVectorMap populations2;
        processing::getPopulations( populations2, points );
        PidT sum(0);
        for ( auto it = populations2.begin(); it != populations2.end(); ++it )
        {
            if ( it->first == PrimitiveT::LONG_VALUES::UNSET )
                continue;

            if ( it->first < 10 )
                std::cout << "it:" << it->first << ", " << it->second.size() << ", sum: " << sum;
            sum+= it->second.size();
            if ( it->first < 10 )
                std::cout << " -> " << sum << std::endl;
        }
        std::cout << "reduction: " << points.size() << " --> " << sum << " == " << Scalar(sum) / points.size() << std::endl;

        for ( auto it = gidsBySize.begin(); it != gidsBySize.end(); ++it )
        {
            const PidT population = it->first;
            const GidT gid = it->second;
            if ( population < popLimit ) continue;
            std::cout << "oldPop[" << gid << "]: " << populations[gid].size() << ", new: " << populations2[gid].size() << std::endl;
        }


        std::string input_prims_path = parsePrimitivesPath(argc, argv);
        std::string outName = boost::filesystem::path( input_prims_path).stem().string();
        std::stringstream ssPrims; ssPrims << outName << ".sub_" << ratio << "_" << primLimit << ".csv";
        GF2::io::savePrimitives<PrimitiveT,typename InnerContainerT::const_iterator>( out, ssPrims.str() );

        std::stringstream ssAssoc; ssAssoc << "points_" << outName << ".sub_" << ratio << "_" << primLimit << ".csv";
        GF2::io::writeAssociations<PointPrimitiveT>( points, ssAssoc.str() );

        std::cout << "results written to " << ssPrims.str() << " and " << ssAssoc.str() << "\n";
        std::cout << "../runGlobfit.py -s " << params.scale << " -p " << ssPrims.str() << " -a " << ssAssoc.str() << "\n";

        return EXIT_SUCCESS;
    } //...subsamplePrimitives()
} //...ns globopt

#endif // GO_SUBSAMPLEPRIMITIVES_HPP
