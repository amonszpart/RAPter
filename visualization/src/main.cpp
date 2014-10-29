#include <iostream>
#include "pcl/console/parse.h"

#include "globfit2/visualization/visualization.h"
#include "globfit2/globOpt_types.h"

int main(int argc, char *argv[])
{
    if ( pcl::console::find_switch(argc,argv,"--show") )
    {
        return GF2::vis::showCli<GF2::_2d::PrimitiveT>( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--show3D") )
    {
        return GF2::vis::showCli<GF2::_3d::PrimitiveT>( argc, argv );
    }
    else
    {
        std::cerr << "[" << __func__ << "]: " << "unrecognized cli option" << std::endl;
        std::cout << "usage: " << "--show[3D] --help" << std::endl;
    }
}
