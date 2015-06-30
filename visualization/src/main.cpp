#include <iostream>
#include "pcl/console/parse.h"

#include "rapter/typedefs.h"
#include "rapter/visualization/visualization.h"

#include "rapter/visualization/mst.hpp"
#include "rapter/primitives/impl/planePrimitive.hpp"
#include "rapter/primitives/impl/linePrimitive.hpp"

int main(int argc, char *argv[])
{
    //return rapter::mstMain<float>();

    if ( pcl::console::find_switch(argc,argv,"--show") )
    {
        return rapter::vis::showCli<rapter::_2d::PrimitiveT>( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--show3D") )
    {
        return rapter::vis::showCli<rapter::_3d::PrimitiveT>( argc, argv );
    }
    else
    {
        std::cerr << "[" << __func__ << "]: " << "unrecognized cli option" << std::endl;
        std::cout << "usage: " << "--show[3D] --help" << std::endl;
    }
}
