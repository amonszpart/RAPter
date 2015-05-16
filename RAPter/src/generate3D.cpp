#include "rapter/typedefs.h"                        // _2d, _3d namespaces
#include "rapter/util/parse.h"                      // rapter::console

#include "rapter/optimization/candidateGenerator.h"
#include "rapter/primitives/impl/planePrimitive.hpp"

int generate3D( int argc, char** argv )
{
    if ( rapter::console::find_switch(argc,argv,"--generate3D") )
    {
        return rapter::CandidateGenerator::generateCli< rapter::_3d::PrimitiveContainerT
                                                   , rapter::PointContainerT
                                                   , rapter::Scalar
                                                   , rapter::PointPrimitiveT
                                                   , rapter::_3d::PrimitiveT
                                                   >( argc, argv );
    }
    else
        std::cerr << "[" << __func__ << "]: " << "switched to wrong place..." << std::endl;

    return EXIT_FAILURE;
}
