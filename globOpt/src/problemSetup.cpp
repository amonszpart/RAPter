#include "globfit2/globOpt_types.h" // _2d, _3d namespaces
#include "globfit2/util/parse.h" // GF2::console

#include "globfit2/optimization/problemSetup.h"

int formulate( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--formulate3D") )
    {
        return GF2::ProblemSetup::formulateCli< GF2::_3d::PrimitiveContainerT
                                              , GF2::PointContainerT
                                              , GF2::_3d::PrimitiveT
                                              , GF2::PointPrimitiveT
                                              >( argc, argv );
    }
    else
    {
        return GF2::ProblemSetup::formulateCli< GF2::_2d::PrimitiveContainerT
                                              , GF2::PointContainerT
                                              , GF2::_2d::PrimitiveT
                                              , GF2::PointPrimitiveT
                                              >( argc, argv );
    } //...if find_switch
}
