#include <iostream>
#include "globfit2/globOpt_types.h"
#include "globfit2/io/ps.hpp"
#include "globfit2/io/io.h"
#include "globfit2/util/parse.h"

namespace GF2
{

    template <class _PrimitiveContainerT, class _PointContainerT>
    inline int drawGraphs( int argc, char **argv )
    {
        typedef typename _PrimitiveContainerT::value_type InnerContainerT;
        typedef typename InnerContainerT::value_type PrimitiveT;
        typedef typename _PointContainerT::value_type PointPrimitiveT;

        bool valid = true;
        std::string primPath, cloudPath = "./cloud.ply", assPath, outPath = "graphNew.gv";

        if ( GF2::console::parse_argument( argc, argv, "-p", primPath) < 0 )
            valid = false;
        if ( GF2::console::parse_argument( argc, argv, "-a", assPath) < 0 )
            valid = false;
        if ( GF2::console::parse_argument( argc, argv, "-c", cloudPath) < 0
             && !boost::filesystem::exists(cloudPath) )
            valid = false;
        GF2::console::parse_argument( argc, argv, "-o", outPath );

        if ( !valid )
        {
            std::cout << "Usage: " << argv[0] << " -p prims.csv -a points_prims.csv -c cloud.ply [-o out.gv] [-s scale] [--dids did0,did1]" << std::endl;
            return EXIT_FAILURE;
        }

        _PrimitiveContainerT primitives;
        containers::PrimitiveContainer<PrimitiveT> primsMap;
        io::readPrimitives<PrimitiveT,InnerContainerT>( primitives, primPath, &primsMap );

        _PointContainerT points;
        GF2::io::readPoints<PointPrimitiveT>( points, cloudPath, NULL );

        // read assoc
        {
            std::vector< std::pair<PidT,LidT> > points_primitives;
            GF2::io::readAssociations( points_primitives, assPath, NULL );
            for ( size_t i = 0; i != points.size(); ++i )
            {
                // store association in point
                points[i].setTag( PointPrimitiveT::TAGS::GID, points_primitives[i].first );
            }
        }

        bool showClusters = GF2::console::find_switch( argc, argv, "--clusters" );
        bool writeNames   = GF2::console::find_switch( argc, argv, "--names" );

        std::vector<DidT> dids;
        GF2::console::parse_x_arguments( argc, argv, "--dids", dids );
        std::vector< GidT > tmpEdgeSources;
        GF2::console::parse_x_arguments( argc, argv, "--edge-sources", tmpEdgeSources );
        std::vector< std::pair<GidT,DidT> > edgeSources;
        for ( int i = 0; i < tmpEdgeSources.size(); i += 2 )
            edgeSources.push_back( std::pair<GidT,DidT>(tmpEdgeSources[i], tmpEdgeSources[i+1] ) );


        Scalar pw( 0. );
        GF2::console::parse_argument( argc, argv, "--pw", pw );
        Scalar radius = .1;
        GF2::console::parse_argument( argc, argv, "--radius", radius ); // .1, .3
        Scalar subSample = 1.;
        GF2::console::parse_argument( argc, argv, "--subsample", subSample ); // 0..1

        bool drawOld = GF2::console::find_switch( argc, argv, "--old" );
        bool colourCloud = GF2::console::find_switch( argc, argv, "--colourCloud" );

        io::drawGraph( primsMap, points, outPath, drawOld, true, dids.size() ? &dids : NULL, pw, showClusters, edgeSources.size() ? &edgeSources : NULL );

        Scalar scale = 0.02f;
        if ( GF2::console::parse_argument( argc, argv, "-s", scale ) < 0 )
            std::cerr << "[" << __func__ << "]: " << "can't do drawPs, need -s scale" << std::endl;
        else
        {
            io::drawPs( primsMap, points, outPath + ".ps", scale
                      , /* show: */ true
                      , /* writeNAmes: */ writeNames
                      , colourCloud
                      , &dids
                      , subSample
                      , radius
                      );
        }
    }
}

int main( int argc, char *argv[] )
{
    bool do3D = GF2::console::find_switch( argc, argv, "--3D" );

    if ( do3D )
    {
        return GF2::drawGraphs< GF2::_3d::PrimitiveContainerT
                , GF2::PointContainerT
                >( argc, argv );
    }
    else
        return GF2::drawGraphs< GF2::_2d::PrimitiveContainerT
                , GF2::PointContainerT
                >( argc, argv );

    return 0;
} //...main
