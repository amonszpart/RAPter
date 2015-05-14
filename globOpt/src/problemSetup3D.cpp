#include "rapter/globOpt_types.h" // _2d, _3d namespaces
#include "rapter/util/parse.h" // rapter::console

#include "rapter/optimization/problemSetup.h"
#include "rapter/optimization/impl/problemSetup.hpp"
#include "rapter/primitives/impl/planePrimitive.hpp"

int formulate3D( int argc, char** argv )
{
    if ( rapter::console::find_switch(argc,argv,"--formulate3D") )
    {
        return rapter::ProblemSetup::formulateCli< rapter::_3d::PrimitiveContainerT
                                              , rapter::PointContainerT
                                              , rapter::_3d::PrimitiveT
                                              , rapter::PointPrimitiveT
                                              , rapter::_3d::MyFinitePlaneToFinitePlaneCompatFunctor
                                              >( argc, argv );
    }
    else
    {
        std::cerr << "[" << __func__ << "]: " << "switch error" << std::endl;
        return EXIT_FAILURE;
    }
}

