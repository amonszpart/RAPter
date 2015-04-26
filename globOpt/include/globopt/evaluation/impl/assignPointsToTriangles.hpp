#ifndef GO_ASSIGNPOINTSTOTRIANGLES_HPP
#define GO_ASSIGNPOINTSTOTRIANGLES_HPP

#include "globfit2/simple_types.h"
#include "globfit2/io/inputParser.hpp"
#include "globopt/io/trianglesFromObj.h"
#include "globopt/primitives/impl/triangle.hpp"
#include "pcl/PolygonMesh.h"

#include "pcl/visualization/pcl_visualizer.h"
#include <chrono>
#include <iomanip> // write to filpe

namespace globopt
{

    template <class _PrimitiveContainerT, class _PointContainerT>
    inline std::string writePrimAssocCloud( _PrimitiveContainerT const& prims, _PointContainerT const& points, std::string stem, std::string dir = "./" )
    {
        typedef typename _PointContainerT::PrimitiveT PointPrimitiveT;

        char outPrimsPath[ 512 ], outAssocPath[512], outCloudPath[512];
        sprintf( outPrimsPath, "%s/%s.csv", dir.c_str(), stem.c_str() );
        std::cout << " writing " << prims.size() << " primitives to " << outPrimsPath << std::endl;
        GF2::io::savePrimitives<typename _PrimitiveContainerT::PrimitiveT,typename _PrimitiveContainerT::mapped_type::const_iterator >( prims, std::string(outPrimsPath) );

        sprintf( outAssocPath, "%s/points_%s.csv", dir.c_str(), stem.c_str() );
        std::cout << " writing " << points.size() << " assignments to " << outAssocPath << std::endl;
        GF2::io::writeAssociations<PointPrimitiveT>( points, outAssocPath );

        sprintf( outCloudPath, "%s/%s.cloud.ply", dir.c_str(), stem.c_str() );
        std::cout << " writing " << points.size() << " points to " << outCloudPath << std::endl;
        GF2::io::writePoints<PointPrimitiveT>( points, outCloudPath );

        std::stringstream ss;
        ss << "../show.py -p " << outPrimsPath << " -a " << outAssocPath << " --cloud " << outCloudPath;

        return ss.str();
    }

    template <typename _PointContainerT, typename _PrimitivesMapT, class _PopulationsT, typename _Scalar>
    inline GF2::GidT getClosestPrimitive( _PointContainerT const& points, GF2::PidT const pId, _PrimitivesMapT & primitives, _PopulationsT const& populations, _Scalar scale )
    {
        typedef typename _PrimitivesMapT::PrimitiveT PrimitiveT;
        typedef typename _PointContainerT::PrimitiveT PointPrimitiveT;

        PointPrimitiveT const& point = points[ pId ];

        _Scalar    min_dist( std::numeric_limits<_Scalar>::max() ), tmp;
        GF2::GidT min_gid ( PrimitiveT::LONG_VALUES::UNSET );
        typename PrimitiveT::ExtremaT extrema;
        for ( typename _PrimitivesMapT::ConstIterator it(primitives); it.hasNext(); it.step() )
        {
            auto popIt = populations.find(it.getGid());
            if ( popIt == populations.end() )
                continue;
#           pragma omp critical (CRIT_GETEXTENT)
            {
                it->template getExtent<PointPrimitiveT>( extrema, points, scale, &(popIt->second), /* force_axis_aligned: */ true );
            }
            tmp = it->getFiniteDistance( extrema, point.template pos() );

            if ( tmp < min_dist )
            {
                min_dist = tmp;
                min_gid = it.getGid();
            }
        }
        return min_gid;
    } //...closestPrimitive

    template <typename TriangleT>
    inline void addTriangle( pcl::visualization::PCLVisualizer::Ptr &vptr, TriangleT triangle, GF2::LidT id, Eigen::Vector3f colour = Eigen::Vector3f::Ones() )
    {
        std::stringstream triangleName;
        triangleName << "triangle" << id;
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr polyCloud( new pcl::PointCloud<pcl::PointXYZRGB> );
        for ( int j = 0; j != 3; ++j )
        {
            pcl::PointXYZRGB pnt; pnt.getVector3fMap() = triangle.getCorner(j);
            polyCloud->push_back( pnt );
        }

        vptr->removePolygonMesh( triangleName.str() );
        vptr->addPolygon<pcl::PointXYZRGB>( polyCloud, colour(0), colour(1), colour(2), triangleName.str(), 0 );
    }

    // Release/bin/eval --planes --assign gt.obj --cloud cloud.ply --scale 0.02 -p primitives.globfit.csv -a points_primitives.globfit.csv --silent --filter-ambig 1
    template < class _PrimitiveVectorT
             , class _PrimitiveMapT
             , class _PointContainerT
             , class _PclCloudT>
    int assignPointsToTriangles( int argc, char** argv )
    {
        using namespace GF2;

        typedef typename _PrimitiveMapT::mapped_type    InnerContainerT;
        typedef typename _PclCloudT::Ptr                PclPtrT;
        typedef typename InnerContainerT::value_type    PrimitiveT;
        typedef typename PrimitiveT::Scalar             Scalar;
        typedef typename _PointContainerT::PrimitiveT   PointPrimitiveT;
        typedef typename _PointContainerT::const_iterator   PointConstIteratorT;
        typedef typename Eigen::Matrix<Scalar,3,1> Vector;
        typedef Triangle<Scalar> Triangle;

        _PointContainerT    points;
        PclPtrT             pclCloud;
        _PrimitiveVectorT   primitivesVector;
        _PrimitiveMapT      primitives;
        struct MyParams { Scalar scale; } params;

        std::string meshPath("");
        GF2::console::parse_argument( argc, argv, "--assign", meshPath );
        if ( !boost::filesystem::exists(meshPath) )
        {
            std::cerr << "[" << __func__ << "]: " << "mesh at \"" << meshPath <<  "\" does not exist, please specify with --assign x.ply/obj" << std::endl;
            return EXIT_FAILURE;
        }

        // read input
        GF2::parseInput<InnerContainerT,_PclCloudT>( points, pclCloud, primitivesVector, primitives, params, argc, argv, /* readAssoc: */ true );

        bool silent = GF2::console::find_switch( argc, argv, "--silent" );

        int ambig(0);
        GF2::console::parse_argument( argc, argv, "--filter-ambig", ambig );
        std::cout << "[" << __func__ << "]: " << "Filtering points, that have more than " << ambig << " triangles inside " << params.scale << " radius\n";

        Scalar minPlaneEdge( 0 );
        GF2::console::parse_argument( argc, argv, "--min-plane-edge", minPlaneEdge );
        std::cout << "[" << __func__ << "]: " << "Filtering triangles and points, that have an edge shorter, than " << minPlaneEdge << "\n";

        LidT N(10000);
        GF2::console::parse_argument( argc, argv, "--n-rels", N );
        std::cout << "[" << __func__ << "]: " << "Using " << N << " pairwise comparisons, change with --n-rels" << std::endl;

        Scalar recallPlaneDistThreshold(params.scale);
        GF2::console::parse_argument( argc, argv, "--recall-thresh", recallPlaneDistThreshold );
        std::cout << "[" << __func__ << "]: " << "Valid assignment, if plane closer to GT than " << recallPlaneDistThreshold << ", change with --recall-thresh" << std::endl;

        std::string statLogPath( "" );
        GF2::console::parse_argument( argc, argv, "--stat-log", statLogPath );
        std::cout << "[" << __func__ << "]: " << "Appending stats to \"" << statLogPath << "\", change with --stat-log" << std::endl;

        // parse input mesh
        std::vector<Triangle> triangles;
        getTrianglesFromObj( triangles, meshPath, minPlaneEdge );
        std::cout << "have " << triangles.size() << " triangles" << std::endl;

        // debug ( add triangles )
//        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
//        vptr->setBackgroundColor( .5, .6, .6 );
//        for ( int i = 0; i != triangles.size(); ++i )
//        {
//            addTriangle( vptr, triangles[i], i );
//        }

//        // debug
//        vptr->addPointCloud<typename _PclCloudT::PointType>( pclCloud );
//        vptr->addCoordinateSystem( 0.1 );
//        vptr->spinOnce(200);
//        vptr->spin();

        // parse input names
        std::string primPath        = parsePrimitivesPath( argc, argv );
        std::cout << "primPath.rfind(/): " << primPath.rfind("/") << std::endl;
        std::string primStem        = primPath;
        if ( primPath.find("/") != std::string::npos )
            primStem = primPath.substr( primPath.rfind("/")+1 );
        std::cout << "primStem0: " << primStem << std::endl;
        primStem = primStem.substr( 0, primStem.size() - 3 );
        std::cout << "primStem1: " << primStem << std::endl;
        std::string outAnglesPath   = primStem + "angles.csv";
        std::string outAnglesSimplePath = primStem + "simple.angles.csv";
        std::string assocPath       = parseAssocPath( argc, argv );
        std::string assocStem       = assocPath.substr( 0, assocPath.size() - 3 );

        // work: 1. match point to triangles                         --> [points_]primitives.gt.csv
        //       2. compare triangles to already assigned primitives -->                    .angles.csv

        GF2::GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        GF2::processing::getPopulations( populations, points );

        // pointid => < triangleId, primitiveGid >
        PidT reassignedCount = 0;
        std::map<PidT,std::pair< LidT, LidT> > pointsTriangles;
        PidT unambigGtPointsCount = 0; // number of points, that have a triangle assigned
        // pointid
        //PidT pId( 0 );
#       pragma omp parallel for num_threads(6)
        for ( PidT pId = 0; pId < points.size(); ++pId )
        {
            // cache point reference
            //PointPrimitiveT const& point = *pIt;
            PointPrimitiveT const& point = points[ pId ];
            // skip not assigned points
            //if ( point.gidUnset() ) continue;

            // cache primitive
            //PrimitiveT const& prim = primIt->second.at(0);

            // select closest triangle
            Scalar  minPointTriangleDistance( std::numeric_limits<Scalar>::max() );
            GidT    closestTriangleId       ( -1 ),
                    triangleId              ( 0  );
            LidT    trianglesNearby         ( 0 );
            // cache point position
            Vector  pos                     ( point.template pos() );
            Vector  triangleNormal;
            // iterate triangles
            for ( auto triIt = triangles.begin(); triIt != triangles.end(); ++triIt, ++triangleId )
            {
                Scalar dist = triIt->getDistance( pos );
                // note, if closer and close enough
                if ( (dist < minPointTriangleDistance) && (dist < params.scale) )
                {
                    minPointTriangleDistance = dist;
                    closestTriangleId        = triangleId;
                    Scalar triangleNormalAngle = 0.;
                    if ( !trianglesNearby )
                        triangleNormal = triIt->dir();
                    else
                    {
                        triangleNormalAngle = GF2::angleInRad( triangleNormal, triIt->dir() );
                        triangleNormalAngle = std::min( triangleNormalAngle, Scalar(M_PI) - triangleNormalAngle );
                    }

                    ++trianglesNearby;

                    if ( ambig && (trianglesNearby > ambig) && (triangleNormalAngle > 0.0001) )
                    {
                        closestTriangleId = -1;
                        break;
                    } //...if ambiguousity threshold exceeded
                }
            } //...for triangles

            // if triangle found
            if ( (closestTriangleId >= 0) )
            {
#               pragma omp critical (UNAMBCNT)
                {
                    ++unambigGtPointsCount;
                }

                // we *need* a primitive for this point, since it ended up in the GT
                GidT gid( PrimitiveT::LONG_VALUES::UNSET );

                // get point group Id
                gid = point.getTag( PointPrimitiveT::TAGS::GID );

                // make sure primitive exists
                if ( gid != PrimitiveT::LONG_VALUES::UNSET )
                {
                    auto primIt = primitives.find( gid );
                    if ( primIt == primitives.end() || !primIt->second.size() )
                        gid = PrimitiveT::LONG_VALUES::UNSET;
                }

                // assign closest, if not assigned
                if ( gid == PrimitiveT::LONG_VALUES::UNSET )
                {
                    gid = getClosestPrimitive( points, pId, primitives, populations, params.scale );
                    ++reassignedCount;
                }

                // remember point for later
                if ( gid != PrimitiveT::LONG_VALUES::UNSET )
                {
#                   pragma omp critical (PTR)
                    {
                        // note point to triangle assignment
                        pointsTriangles[ pId ] = std::pair<LidT,GidT>( closestTriangleId, gid );
                    }
                } //...gid not unset

            } //...triangle found
            else if ( !silent )
            {
                std::cerr << "no triangle found for point..." << std::endl;
            }

            if ( !(pId % 100000) )
                std::cout << pId << std::endl;
        } //...points
        std::cout << "[" << __func__ << "]: " << "reassigned " << reassignedCount << " points" << std::endl;

        // store new assignments
        _PointContainerT gtPoints;
        gtPoints.reserve( points.size() );
        for ( PidT pid = 0; pid != points.size(); ++pid )
        {
            gtPoints.push_back( PointPrimitiveT(points[pid].template pos(), points[pid].template dir()) );
            gtPoints.back().setTag( PointPrimitiveT::TAGS::PID, pid );
            // clear previous assignment
            //gtPoints.back().setTag( PointPrimitiveT::TAGS::GID, PointPrimitiveT::LONG_VALUES::UNSET );
        } // ...for all points

        // save triangles as primitives
        _PrimitiveMapT gtPrims;
        for ( int triangleId = 0; triangleId != triangles.size(); ++triangleId )
        {
            // plane from triangle centroid
            gtPrims[triangleId].push_back( PrimitiveT(triangles[triangleId].getMean(), triangles[triangleId].template dir()) );
            // assign triangle id
            gtPrims[triangleId].back().setTag( PrimitiveT::TAGS::GID    , triangleId );
            gtPrims[triangleId].back().setTag( PrimitiveT::TAGS::DIR_GID, triangleId );
        }

        std::ofstream fAnglesSimple( outAnglesSimplePath );
        // write angles to file
        _PointContainerT orientedPoints, orientedGtPoints;
        orientedPoints  .reserve( pointsTriangles.size() );
        orientedGtPoints.reserve( pointsTriangles.size() );
        for ( auto it = pointsTriangles.begin(); it != pointsTriangles.end(); ++it )
        {
            // read
            PidT                   pid        = it->first;
            LidT                   triangleId = it->second.first;
            GidT                   primGid    = it->second.second;
            Triangle        const& triangle   = triangles .at( triangleId );
            PrimitiveT      const& prim       = primitives.at( primGid ).at( 0 );
            PointPrimitiveT const& point      = points    .at( pid );

            if ( gtPrims.find(triangleId) == gtPrims.end() )
            {
                std::cerr << "nonexistent triangle..." << std::endl;
                continue;
            }
            PrimitiveT const& gtPrim = gtPrims.at(triangleId).at(0);
#if 1
            if ( (gtPrim.projectPoint(point.template pos()) - prim.projectPoint(point.template pos())).norm() > recallPlaneDistThreshold )
                continue;
#else
            // compare infinite distance to origin
            double dist = std::abs(gtPrims.at(triangleId).at(0).getDistance(Vector::Zero()) - prim.getDistance( Vector::Zero()));
            if ( dist > recallPlaneDistThreshold )
            {
                // 10,23,16,34
                if ( (primGid == 10) || (primGid == 16) || (primGid == 23) || (primGid == 34) )
                {
                    std::cout << "dist: " << dist << std::endl;
                    std::cout << "prim[" << primGid << "]: " << prim.toString() << std::endl;
                    std::cout << "gt: " << gtPrims.at(triangleId).at(0).toString() << std::endl;
                    std::cout << "triangle: " << triangles[triangleId].template dir().transpose() <<std::endl;
                    pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
                    vptr->setBackgroundColor( .5, .6, .6 );
                    PrimitiveT::template draw<PointPrimitiveT>( /*   primitive: */ gtPrims.at(triangleId).at(0)
                                                              , /*      points: */ points
                                                              , /*   threshold: */ params.scale
                                                              , /*     indices: */ &populations.at(primGid)
                                                              , /*      viewer: */ vptr
                                                              , /*   unique_id: */ "gtPrim"
                                                              , /*      colour: */ 0,1,0
                                                              , /* viewport_id: */ 0
                                                              , /*     stretch: */ 1.0
                                                              , /*       qhull: */ 0
                                                              );
                    //vptr->addArrow( b, a, 0,1,0, true, "gtPrim",)
                    PrimitiveT::template draw<PointPrimitiveT>( /*   primitive: */ prim
                                                              , /*      points: */ points
                                                              , /*   threshold: */ params.scale
                                                              , /*     indices: */ &populations.at(primGid)
                                                              , /*      viewer: */ vptr
                                                              , /*   unique_id: */ "prim"
                                                              , /*      colour: */ 1,0,0
                                                              , /* viewport_id: */ 0
                                                              , /*     stretch: */ 1.0
                                                              , /*       qhull: */ 0
                                                              );
                    vptr->spin();
                }
                continue;
            }
#endif

            orientedPoints.push_back( PointPrimitiveT(point.template pos(),prim.template dir()) );
            orientedPoints.back().setTag( PointPrimitiveT::TAGS::GID, primGid );
            orientedPoints.back().setTag( PointPrimitiveT::TAGS::PID, orientedPoints.size()-1 );

            orientedGtPoints.push_back( PointPrimitive(point.template pos(),triangle.template dir()) );

            // save point assignment
            gtPoints.at( pid ).setTag( PointPrimitiveT::TAGS::GID, triangleId );

            // Deprecated
            // calculate angle
            Scalar angle = GF2::angleInRad( prim.template dir(), triangle.template dir() );
            // reduce to parallel
            angle = std::min( angle, Scalar(std::abs(M_PI-angle)));

            // write angle
            fAnglesSimple << std::setprecision(9) << angle << "\n";
        } //...for all points that have triangles
        fAnglesSimple.close();
#if 1
        double avg = 0.;
        PidT   pId0, pId1;
        //double angle(0.);
        std::ofstream fAngles( outAnglesPath );

        auto angleLambda = [&orientedPoints, &orientedGtPoints, &avg, &fAngles, &N]( PidT const pId0, PidT const pId1 )
        {
            double angle = std::abs( GF2::angleInRad( orientedPoints  [pId0].template dir().template cast<double>(), orientedPoints  [pId1].template dir().template cast<double>() )
                                   - GF2::angleInRad( orientedGtPoints[pId0].template dir().template cast<double>(), orientedGtPoints[pId1].template dir().template cast<double>() ) );
            while ( angle > M_PI ) { angle -= M_PI; std::cout << "hit" << std::endl; }
            angle = std::min( angle, double(M_PI - angle) );
            //outAngles.push_back( angle );
            avg += angle / double(N);
            fAngles << std::setprecision(9) << angle << "\n";
        };

        if ( N > orientedPoints.size() * orientedPoints.size() / 2. )
        {
            std::cout << "[" << __func__ << "]: " << "exhaustive" << std::endl;
            for ( PidT pId0 = 0; pId0 != orientedPoints.size()-1; ++pId0 )
                for ( PidT pId1 = pId0+1; pId1 != orientedPoints.size(); ++pId1 )
                    angleLambda( pId0, pId1 );
        }
        else
        {
            std::cout << "[" << __func__ << "]: " << "randomized: " << N << " < " << orientedPoints.size() * orientedPoints.size() * 0.5 << "(" << points.size() << " * " << points.size() << "=" << points.size() * points.size() << std::endl;
            for ( PidT id0 = 0; id0 != N; ++id0 )
            {
                do {
                    pId0 = rand() % orientedPoints.size();
                    pId1 = rand() % orientedPoints.size();
                }
                while (pId0 == pId1);

                angleLambda( pId0, pId1 );
            }
        }

        std::cout << "\n\n mean: " << avg << " rad, " << avg * 180.0 / M_PI << " deg" << std::endl;
        double recall = orientedPoints.size() / double(unambigGtPointsCount); //Scalar(pointsTriangles.size());
        std::cout << "recall: " << recall * Scalar(100.0) << " %" << std::endl << "\n\n";
        if ( !statLogPath.empty() )
        {
            std::ofstream fStats( statLogPath, std::ios_base::app );
            fStats << recall * 100.0 << "," << avg << std::endl;
            fStats.close();
            std::cout << "appended recall,avg to " << statLogPath << std::endl;
        }

        // close file
        fAngles.close();
#else
        // subsample
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::vector<LidT> ids( orientedPoints.size() );
        {
            for ( LidT id = 0; id != orientedPoints.size(); ++id )
                ids[id] = id;
            shuffle( ids.begin(), ids.end(), std::default_random_engine(seed) );
        }

        // pairwise relative comparisons
        Scalar angle;
        std::cout << "have " << orientedPoints.size() << " orientedPOints and " << orientedGtPoints.size() << " gt points\n" << std::endl;
        if ( orientedPoints.size() * orientedPoints.size() / 2. < N )
        {
            N = 0.5 * orientedPoints.size() * orientedPoints.size();
            std::cout << " set N to " << N << std::endl;
        }
        //for ( PidT pId0 = 0; pId0 != orientedPoints.size() - 1; ++pId0 )
        //    for ( PidT pId1 = pId0+1; pId1 != orientedPoints.size(); ++pId1 )
        PidT sqrtN = std::floor( std::sqrt(Scalar(2.*N)) );
        std::cout << "sqrtN: " << sqrtN << ", N*N-1/2: " << sqrtN * (sqrtN-1) * 0.5 << ", N: " << N << std::endl;
        //std::vector<Scalar> outAngles; outAngles.reserve(N);
        PidT pId0, pId1;
        //Scalar mean = 0.;
        //size_t meanCount = 0;
        std::ofstream fAngles( outAnglesPath );
        for ( PidT id0 = 0; id0 != sqrtN-1; ++id0 )
        {
            if ( id0 % 10000 == 0 )
                std::cout << "id0 " << id0 << " / " << sqrtN-1 << std::endl;
            for ( PidT id1 = id0+1; id1 != sqrtN; ++id1/*, ++meanCount*/ )
            {
                pId0 = ids[id0];
                pId1 = ids[id1];
                angle = std::abs( GF2::angleInRad( orientedPoints  [pId0].template dir(), orientedPoints  [pId1].template dir() )
                                - GF2::angleInRad( orientedGtPoints[pId0].template dir(), orientedGtPoints[pId1].template dir() ) );
                while ( angle > M_PI ) angle -= M_PI;
                angle = std::min( angle, Scalar(M_PI - angle) );
                //outAngles.push_back( angle );
                //mean += angle;
                fAngles << angle << "\n";
            }
        }

        // close file

        fAngles.close();
#endif
        // log
        std::cout << "[" << __func__ << "]: " << "wrote angles to " << outAnglesPath << std::endl;
        std::cout << "written to " << outAnglesSimplePath << std::endl;

        // write primitives
        std::stringstream ssPrims; ssPrims << primStem << "gt.csv";
        GF2::io::savePrimitives<PrimitiveT, typename InnerContainerT::const_iterator>( gtPrims, ssPrims.str() );

        // write association
        std::stringstream ssAssoc; ssAssoc << assocStem << "gt.csv";
        GF2::io::writeAssociations<PointPrimitiveT>( gtPoints, ssAssoc.str() );
        std::cout << "results written to " << ssPrims.str() << " and " << ssAssoc.str() << "\n";
        std::cout << "../show.py --cloud -s " << params.scale << " -p " << ssPrims.str() << " -a " << ssAssoc.str() << " --save-poly" << std::endl;

        std::stringstream cmd;
        cmd << writePrimAssocCloud( primitives, orientedPoints, primStem + "recalled" );
        cmd << " -s " << params.scale;
        std::cout << cmd.str() << std::endl;

        int err = system( (cmd.str() + " &").c_str() );
        if ( err ) std::cout << "err: " << err << std::endl;

        return EXIT_SUCCESS;
    } //...assignPointsToTriangles

} //...ns globopt

#endif // GO_ASSIGNPOINTSTOTRIANGLES_HPP

