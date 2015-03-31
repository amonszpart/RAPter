#include "globfit2/comparison/globFit.h"

#include <iostream>
#include "globfit2/util/parse.h"
#include "globfit2/globOpt_types.h"

namespace globopt
{
    namespace globfit
    {
        inline int printUsage( int argc, char** argv )
        {
            std::cout << "usage: --planes or --lines " << std::endl;
            return EXIT_SUCCESS;
        }
    }
}

int main(int argc, char *argv[])
{
    if ( GF2::console::find_switch(argc,argv,"--help") )
    {
        return globopt::globfit::printUsage( argc, argv );
    }
    else if ( GF2::console::find_switch(argc,argv,"--planes") )
    {
        if ( GF2::console::find_switch(argc,argv,"--from") )
            return globopt::fromGlobFit< GF2::PclCloudT
                                       , GF2::_3d::PrimitiveMapT
                                       , GF2::_3d::PrimitiveVectorT
                                       , GF2::PointContainerT
                                       >( argc, argv );
        else
            return globopt::toGlobFit< GF2::PclCloudT
                                     , GF2::_3d::PrimitiveMapT
                                     , GF2::_3d::PrimitiveVectorT
                                     , GF2::PointContainerT
                                     >( argc, argv );
    }
    else if ( GF2::console::find_switch(argc,argv,"--lines") )
    {
        std::cout << "--lines (2D) unimplemented" << std::endl;
        return EXIT_FAILURE;
    }
    else
        return globopt::globfit::printUsage( argc, argv );

    return EXIT_FAILURE;
}
