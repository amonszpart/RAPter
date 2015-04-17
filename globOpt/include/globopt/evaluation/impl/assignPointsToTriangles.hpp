#ifndef GO_ASSIGNPOINTSTOTRIANGLES_HPP
#define GO_ASSIGNPOINTSTOTRIANGLES_HPP

#include "globfit2/simple_types.h"
#include "globfit2/io/inputParser.hpp"
#include "globopt/io/trianglesFromObj.h"
#include "globopt/primitives/impl/triangle.hpp"
#include "pcl/PolygonMesh.h"

#include "pcl/visualization/pcl_visualizer.h"

namespace globopt
{

    template <typename _Scalar>
    _Scalar sgn( _Scalar a )
    {
        return (a > _Scalar(0.)) - (a < _Scalar(0.));
        //_Scalar sign = (a > _Scalar(0.)) - (a < _Scalar(0.));

        //return (sign == _Scalar(0.)) ? _Scalar(1.) : sign;
    }

#if 0
    template <typename Scalar, typename _VectorT, class _TriangleT>
    inline Scalar getDistance( _TriangleT const& tri, _VectorT const& point )
    {
        _VectorT p1  = tri.getCorner(0);
        _VectorT p2  = tri.getCorner(1);
        _VectorT p3  = tri.getCorner(2);
         // calculate vectors from point f to vertices p1, p2 and p3:
         _VectorT f1 = p1 - point;
         _VectorT f2 = p2 - point;
         _VectorT f3 = p3 - point;

         // calculate the areas (parameters order is essential in this case):
         _VectorT va  = ( p1-p2 ).cross( p1-p3 ); // main triangle cross product
         _VectorT va1 = f2.cross( f3 ); // p1's triangle cross product
         _VectorT va2 = f3.cross( f1 ); // p2's triangle cross product
         _VectorT va3 = f1.cross( f2 ); // p3's triangle cross product
         Scalar   a   = va.norm(); // main triangle area

         // calculate barycentric coordinates with sign:
         Scalar   a1  = va1.norm() / a * sgn( va.dot(va1) );
         Scalar   a2  = va2.norm() / a * sgn( va.dot(va2) );
         Scalar   a3  = va3.norm() / a * sgn( va.dot(va3) );

         // find the uv corresponding to point f (uv1/uv2/uv3 are associated to p1/p2/p3):
         std::cout << "a1: " << a1 << ", a2: " << a2 << ", a3: " << a3 << std::endl;

         return -1.f;
    } // ..getDistance

#endif

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

    // Usage: Release/bin/eval --planes --assign gtHandAligned.obj --cloud cloud.ply -p primitives_it9.bonmin.csv -a points_primitives_it8.csv --scale 0.05
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

        int  ambig(0);
        GF2::console::parse_argument( argc, argv, "--filter-ambig", ambig );
        std::cout << "[" << __func__ << "]: " << "Filtering points, that have more than " << ambig << " triangles inside " << params.scale << " radius\n";

        // parse input mesh
        std::vector<Triangle> triangles;
        getTrianglesFromObj( triangles, meshPath );
        std::cout << "have " << triangles.size() << " triangles" << std::endl;

        // debug ( add triangles )
        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
        vptr->setBackgroundColor( .5, .6, .6 );
        for ( int i = 0; i != triangles.size(); ++i )
        {
            addTriangle( vptr, triangles[i], i );
        }

        // debug
        vptr->addPointCloud<typename _PclCloudT::PointType>( pclCloud );
        vptr->addCoordinateSystem( 0.1 );
        vptr->spinOnce(200);

        // parse input names
        std::string primPath        = parsePrimitivesPath( argc, argv );
        std::string primStem        = primPath.substr( 0, primPath.size() - 3 );
        std::string outAnglesPath   = primStem + "angles.csv";
        std::string assocPath       = parseAssocPath( argc, argv );
        std::string assocStem       = assocPath.substr( 0, assocPath.size() - 3 );

        // work: 1. match point to triangles                         --> [points_]primitives.gt.csv
        //       2. compare triangles to already assigned primitives -->                    .angles.csv

        // pointid => < triangleId, primitiveGid >
        std::map<PidT,std::pair< LidT, LidT> > pointsTriangles;
        // pointid
        PidT pId( 0 );
        for ( PointConstIteratorT pIt = points.begin(); pIt != points.end(); ++pIt, ++pId )
        {
            // cache point reference
            PointPrimitiveT const& point = *pIt;
            // skip not assigned points
            if ( point.gidUnset() ) continue;

            // get point group Id
            const GidT gid = point.getTag( PointPrimitiveT::TAGS::GID );

            // get assigned primitive
            auto primIt = primitives.find( gid );
            // skip, if not found
            if ( primIt == primitives.end() )
            {
                if ( !silent )
                    std::cerr << "could not find primitive by gid " << gid << ", for point " << pId << std::endl; fflush(stderr);
                continue;
            }

            // skip, if patch is empty
            if ( !primIt->second.size() )
            {
                std::cerr << "gid " << gid << " is empty...no primitive in group" << std::endl; fflush(stderr);
                continue;
            }

            // cache primitive
            //PrimitiveT const& prim = primIt->second.at(0);

            // select closest triangle
            Scalar  minPointTriangleDistance( std::numeric_limits<Scalar>::max() );
            GidT    closestTriangleId       ( -1 ),
                    triangleId              ( 0  );
            LidT    trianglesNearby         ( 0 );
            // cache point position
            Vector  pos                     ( point.template pos() );
            // iterate triangles
            for ( auto triIt = triangles.begin(); triIt != triangles.end(); ++triIt, ++triangleId )
            {
                // estimate distance
                Scalar dist = triIt->getDistance( pos );
                // note, if closer and close enough
                if ( (dist < minPointTriangleDistance) && (dist < params.scale) )
                {
                    minPointTriangleDistance = dist;
                    closestTriangleId        = triangleId;
                    ++trianglesNearby;

                    if ( ambig && (trianglesNearby > ambig) )
                    {
                        closestTriangleId = -1;
                        break;
                    } //...if ambiguousity threshold exceeded
                }
            } //...for triangles

            // if triangle found
            if ( (closestTriangleId >= 0) )
            {
                // note point to triangle assignment
                pointsTriangles[ pId ] = std::pair<LidT,GidT>( closestTriangleId, gid );
            }
            else if ( !silent )
            {
                std::cerr << "no triangle found for point..." << std::endl;
            }

            if ( !(pId % 100000) )
                std::cout << pId << std::endl;
        } //...points

        // store new assignments
        _PointContainerT gtPoints;
        gtPoints.reserve( points.size() );
        for ( PidT pid = 0; pid != points.size(); ++pid )
        {
            gtPoints.push_back( points[pid] );
            // clear previous assignment
            gtPoints.back().setTag( PointPrimitiveT::TAGS::GID, PointPrimitiveT::LONG_VALUES::UNSET );
        } // ...for all points

        // write angles to file
        std::ofstream fAngles( outAnglesPath );
        for ( auto it = pointsTriangles.begin(); it != pointsTriangles.end(); ++it )
        {
            // read
            PidT              pid        = it->first;
            LidT              triangleId = it->second.first;
            GidT              primGid    = it->second.second;
            Triangle   const& triangle   = triangles.at( triangleId );
            PrimitiveT const& prim       = primitives.at( primGid ).at( 0 );

            // calculate angle
            Scalar angle = GF2::angleInRad( prim.template dir(), triangle.template dir() );
            // reduce to parallel
            angle = std::min( angle, Scalar(std::abs(M_PI-angle)));

            // write angle
            fAngles << angle << "\n";

            // save point assignment
            gtPoints.at( pid ).setTag( PointPrimitiveT::TAGS::GID, triangleId );
        } //...for all points that have triangles
        // close file
        fAngles.close();
        // log
        std::cout << "[" << __func__ << "]: " << "wrote angles to " << outAnglesPath << std::endl;

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

        // write primitives
        std::stringstream ssPrims; ssPrims << primStem << "gt.csv";
        GF2::io::savePrimitives<PrimitiveT, typename InnerContainerT::const_iterator>( gtPrims, ssPrims.str() );

        // write association
        std::stringstream ssAssoc; ssAssoc << assocStem << "gt.csv";
        GF2::io::writeAssociations<PointPrimitiveT>( gtPoints, ssAssoc.str() );
        std::cout << "results written to " << ssPrims.str() << " and " << ssAssoc.str() << "\n";
        std::cout << "../show.py -s " << params.scale << " -p " << ssPrims.str() << " -a " << ssAssoc.str() << " --save-poly" << std::endl;

        return EXIT_SUCCESS;
    } //...assignPointsToTriangles

} //...ns globopt

#endif // GO_ASSIGNPOINTSTOTRIANGLES_HPP
