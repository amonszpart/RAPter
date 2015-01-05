#include <iostream>

#if 0
    #include "globfit2/globOpt_types.h" // _2d, _3d namespaces
    #include "globfit2/util/parse.h" // GF2::console

    #include "globfit2/optimization/solver.h" // todo: move datafit to datafit.h
#endif

int datafit( int argc, char** argv )
{
#if 0
    if ( GF2::console::find_switch(argc,argv,"--datafit3D") )
    {
        return GF2::Solver::datafit< GF2::_3d::PrimitiveContainerT
                                   , GF2::_3d::InnerPrimitiveContainerT
                                   , GF2::_3d::PrimitiveT
                                   >( argc, argv );
    }
    else
    {
        return GF2::Solver::datafit< GF2::_2d::PrimitiveContainerT
                                   , GF2::_2d::InnerPrimitiveContainerT
                                   , GF2::_2d::PrimitiveT
                                   >( argc, argv );
    } //...if find_switch
#else
    std::cerr << "commented out" << std::endl;
    return EXIT_FAILURE;
#endif
}

