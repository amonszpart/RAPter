#include <iostream>
#include "globfit2/visualization/visualization.h"


int main(int argc, char *argv[])
{
    std::cout << "visok" << std::endl;
    return GF2::vis::lines::showLinesCli( argc, argv );
}
