#include "globfit2/globOpt_types.h" // _2d, _3d namespaces
#include "globfit2/util/parse.h" // GF2::console

#include "globfit2/optimization/solver.h" // todo: move generate to generate.h

int generate( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--generate3D") )
    {
        return GF2::Solver::generateCli< GF2::_3d::PrimitiveContainerT
                                       , GF2::PointContainerT
                                       , GF2::Scalar
                                       , GF2::PointPrimitiveT
                                       , GF2::_3d::PrimitiveT
                                       >( argc, argv );
    }
    else
    {
        return GF2::Solver::generateCli< GF2::_2d::PrimitiveContainerT
                                       , GF2::PointContainerT
                                       , GF2::Scalar
                                       , GF2::PointPrimitiveT
                                       , GF2::_2d::PrimitiveT
                                       >( argc, argv );
    }

    return EXIT_FAILURE;
}
