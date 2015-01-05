#include "globfit2/globOpt_types.h" // _2d, _3d namespaces
#include "globfit2/util/parse.h" // GF2::console

#include "globfit2/optimization/solver.h"

int solve3D( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--solver3D") )
    {
        return GF2::Solver::solve< GF2::_3d::PrimitiveContainerT
                                 , GF2::_3d::InnerPrimitiveContainerT
                                 , GF2::_3d::PrimitiveT
                                 >( argc, argv );
    }
    else
        std::cerr << "[" << __func__ << "]: " << "wrong switch" << std::endl;
}

