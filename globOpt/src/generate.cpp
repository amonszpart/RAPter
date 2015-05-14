#include "rapter/globOpt_types.h" // _2d, _3d namespaces
#include "rapter/util/parse.h" // rapter::console

#include "rapter/optimization/candidateGenerator.h"

int generate( int argc, char** argv )
{
    if ( rapter::console::find_switch(argc,argv,"--generate") )
    {
        return rapter::CandidateGenerator::generateCli< rapter::_2d::PrimitiveContainerT
                                                   , rapter::PointContainerT
                                                   , rapter::Scalar
                                                   , rapter::PointPrimitiveT
                                                   , rapter::_2d::PrimitiveT
                                                   >( argc, argv );
    }
    else
        std::cerr << "[" << __func__ << "]: " << "switched to wrong place " << std::endl;

    return EXIT_FAILURE;
}
