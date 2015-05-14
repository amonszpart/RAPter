#include "rapter/globOpt_types.h" // _2d, _3d namespaces
#include "rapter/util/parse.h" // rapter::console

#include "rapter/optimization/problemSetup.h"
#include "rapter/optimization/impl/problemSetup.hpp"

int formulate( int argc, char** argv )
{
    if ( rapter::console::find_switch(argc,argv,"--formulate") )
    {
        return rapter::ProblemSetup::formulateCli< rapter::_2d::PrimitiveContainerT
                                              , rapter::PointContainerT
                                              , rapter::_2d::PrimitiveT
                                              , rapter::PointPrimitiveT
                                              , rapter::_2d::MyFiniteLineToFiniteLineCompatFunctor
                                              >( argc, argv );
    } //...if find_switch
    else
    {
        std::cerr << "[" << __func__ << "]: " << "switch error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
