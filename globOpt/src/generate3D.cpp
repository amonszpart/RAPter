#include "globfit2/globOpt_types.h" // _2d, _3d namespaces
#include "globfit2/util/parse.h" // GF2::console

#include "globfit2/optimization/candidateGenerator.h"

int generate3D( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--generate3D") )
    {
        return GF2::CandidateGenerator::generateCli< GF2::_3d::PrimitiveContainerT
                                                   , GF2::PointContainerT
                                                   , GF2::Scalar
                                                   , GF2::PointPrimitiveT
                                                   , GF2::_3d::PrimitiveT
                                                   >( argc, argv );
    }
    else
        std::cerr << "[" << __func__ << "]: " << "switched to wrong place..." << std::endl;

    return EXIT_FAILURE;
}
