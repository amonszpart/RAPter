#include "globfit2/globOpt_types.h" // _2d, _3d namespaces
#include "globfit2/util/parse.h" // GF2::console

#include "globfit2/optimization/solver.h"

int solve( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--solver") )
    {
        return GF2::Solver::solve< GF2::_2d::PrimitiveContainerT
                                 , GF2::_2d::InnerPrimitiveContainerT
                                 , GF2::_2d::PrimitiveT
                                 >( argc, argv );
    } //...if find_switch
    else
    {
        std::cerr << "[" << __func__ << "]: " << "wrong switch" << std::endl;
        return EXIT_FAILURE;
    }
}
