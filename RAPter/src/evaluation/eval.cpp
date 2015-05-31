#include <iostream>
#include "rapter/util/parse.h"
#include "rapter/typedefs.h"
#include "rapter/util/pclUtil.h" // pclCloudT
#include "rapter/evaluation/assignPointsToTriangles.h"

namespace rapter
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
    if ( rapter::console::find_switch(argc,argv,"--help") )
    {
        return rapter::eval::printUsage( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--assign") )
    {
        if ( rapter::console::find_switch(argc,argv,"--planes") )
        {
            return rapter::assignPointsToTriangles<rapter::_3d::PrimitiveVectorT
                    , rapter::_3d::PrimitiveMapT
                    , rapter::PointContainerT
                    , rapter::PclCloudT>
                    ( argc, argv );
        }
        else if ( rapter::console::find_switch(argc,argv,"--lines") )
        {
            std::cout << "--lines (2D) unimplemented" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else
        return rapter::eval::printUsage( argc, argv );

    return EXIT_FAILURE;
}
