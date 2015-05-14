#include "rapter/globOpt_types.h" // _2d, _3d namespaces
#include "rapter/util/parse.h" // rapter::console

#include "rapter/optimization/solver.h"

int solve( int argc, char** argv )
{
    if ( rapter::console::find_switch(argc,argv,"--solver") )
    {
        return rapter::Solver::solve< rapter::_2d::PrimitiveContainerT
                                 , rapter::_2d::InnerPrimitiveContainerT
                                 , rapter::_2d::PrimitiveT
                                 >( argc, argv );
    } //...if find_switch
    else
    {
        std::cerr << "[" << __func__ << "]: " << "wrong switch" << std::endl;
        return EXIT_FAILURE;
    }
}
