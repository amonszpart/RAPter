#include <iostream>
#include "globfit2/util/parse.h"
#include "globfit2/globOpt_types.h"
#include "globopt/evaluation/assignPointsToTriangles.h"

namespace globopt
{
namespace eval
{
inline int printUsage( int argc, char** argv )
{
    std::cout << "usage: --planes or --lines AND --assign x.ply/obj" << std::endl;
    std::cout << "angles: Release/bin/eval --planes --assign gt.obj --cloud cloud.ply -p primitives_it9.bonmin.csv -a points_primitives_it8.csv --scale 0.05" << std::endl;
    return EXIT_SUCCESS;
}
}
}

int main(int argc, char *argv[])
{
    if ( GF2::console::find_switch(argc,argv,"--help") )
    {
        return globopt::eval::printUsage( argc, argv );
    }
    else if ( GF2::console::find_switch(argc,argv,"--assign") )
    {
        if ( GF2::console::find_switch(argc,argv,"--planes") )
        {
            return globopt::assignPointsToTriangles<GF2::_3d::PrimitiveVectorT
                    , GF2::_3d::PrimitiveMapT
                    , GF2::PointContainerT
                    , GF2::PclCloudT>
                    ( argc, argv );
        }
        else if ( GF2::console::find_switch(argc,argv,"--lines") )
        {
            std::cout << "--lines (2D) unimplemented" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else
        return globopt::eval::printUsage( argc, argv );

    return EXIT_FAILURE;
}
