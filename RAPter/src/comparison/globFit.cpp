#include "rapter/comparison/globFit.h"

#include <iostream>
#include "rapter/util/parse.h"
#include "rapter/util/pclUtil.h"                    // PclCloudT
#include "rapter/typedefs.h"
#include "rapter/processing/subsamplePrimitives.h"
#include "rapter/processing/assignmentOps.h"

namespace rapter
{
    namespace globfit
    {
        inline int printUsage( int argc, char** argv )
        {
            std::cout << "usage: --planes or --lines or --subsamle-primitives" << std::endl;
            return EXIT_SUCCESS;
        }
    }
}

int main(int argc, char *argv[])
{
    if ( rapter::console::find_switch(argc,argv,"--help") )
    {
        return rapter::globfit::printUsage( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--planes") )
    {
        if ( rapter::console::find_switch(argc,argv,"--from") )
            return rapter::fromGlobFit< rapter::PclCloudT
                                       , rapter::_3d::PrimitiveMapT
                                       , rapter::_3d::PrimitiveVectorT
                                       , rapter::PointContainerT
                                       >( argc, argv );
        else
            return rapter::toGlobFit< rapter::PclCloudT
                                     , rapter::_3d::PrimitiveMapT
                                     , rapter::_3d::PrimitiveVectorT
                                     , rapter::PointContainerT
                                     >( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--lines") )
    {
        std::cout << "--lines (2D) unimplemented" << std::endl;
        return EXIT_FAILURE;
    }
    else if ( rapter::console::find_switch(argc,argv,"--subsample-primitives") )
    {
        return rapter::subsamplePrimitives<rapter::_3d::PrimitiveVectorT
                                           , rapter::_3d::PrimitiveMapT
                                           , rapter::PointContainerT
                                           , rapter::PclCloudT>
                                           ( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--unass-w-planes") )
    {
        return rapter::approxUnassignedWPlanes<rapter::_3d::PrimitiveVectorT
                                           , rapter::_3d::PrimitiveMapT
                                           , rapter::PointContainerT
                                           , rapter::PclCloudT>
                                           ( argc, argv );
    }
    else
        return rapter::globfit::printUsage( argc, argv );

    return EXIT_FAILURE;
}
