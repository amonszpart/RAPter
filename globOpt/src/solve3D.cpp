#include "rapter/globOpt_types.h" // _2d, _3d namespaces
#include "rapter/util/parse.h" // rapter::console

#include "rapter/optimization/solver.h"
#include "rapter/primitives/impl/planePrimitive.hpp"


int solve3D( int argc, char** argv )
{
    if ( rapter::console::find_switch(argc,argv,"--solver3D") )
    {
        return rapter::Solver::solve< rapter::_3d::PrimitiveContainerT
                                 , rapter::_3d::InnerPrimitiveContainerT
                                 , rapter::_3d::PrimitiveT
                                 >( argc, argv );
    }
    else
        std::cerr << "[" << __func__ << "]: " << "wrong switch" << std::endl;

    return EXIT_FAILURE;
}

