#include <iostream>
#include "pcl/console/parse.h"

#include "globfit2/visualization/visualization.h"
#include "globfit2/primitives/linePrimitive2.h"
#include "globfit2/primitives/planePrimitive.h"

int main(int argc, char *argv[])
{
    if ( pcl::console::find_switch("--show") )
    {
        return GF2::vis::showCli<GF2::LinePrimitive2>( argc, argv );
    }
    else if ( pcl::console::find_switch("--show3D") )
    {
        return GF2::vis::showCli<GF2::PlanePrimitive>( argc, argv );
    }
    else
    {
        std::cerr << "[" << __func__ << "]: " << "unrecognized cli option" << std::endl;
    }
}
