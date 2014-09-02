#include <iostream>

#include <pcl/console/parse.h>

#include "optimization/qp/solver.h"

#include "globfit2/optimization/merging.hpp"
#include "globfit2/io/io.h"

int main( int argc, char *argv[] )
{
    if ( (argc == 2) &&
         (   (pcl::console::find_switch(argc,argv,"--help"))
          || (pcl::console::find_switch(argc,argv,"-h"    ))
         )
       )
    {
        std::cout << "[Usage]:\n"
                  << "\t--sample-input\n"
                  << "\t--generate\n"
                  << "\t--formulate\n"
                  << "\t--solver mosek|bonmin|gurobi\n"
                  << "\t--merge\n"
                  << "\t--gfit\n"
                  << "\t--show\n"
                  << std::endl;

        return EXIT_SUCCESS;
    }
    else if ( pcl::console::find_switch(argc,argv,"--sample-input") )
    {
        return GF2::Solver::sampleInput( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--generate") )
    {
        return GF2::Solver::generate( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--formulate") )
    {
        return GF2::Solver::formulate( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--solver") )
    {
        return GF2::Solver::solve( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--gfit") )
    {
        return GF2::Solver::datafit( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--merge") )
    {
        return GF2::Merging::merge<GF2::Solver::PrimitiveContainerT, GF2::Solver::PointContainerT, GF2::Solver::Scalar>( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--show") )
    {
        std::cerr << "[" << __func__ << "]: " << "the show option has been moved to a separate executable, please use thatt one" << std::endl;
        return 1;
        //return GF2::Solver::show( argc, argv );
    }

    return 1;

    // --show --dir . --cloud cloud.ply --scale 0.05f --assoc points_primitives.txt --use-tags --no-clusters --prims primitives.bonmin.txt

//    std::string img_path( "input2.png" );
//    pcl::console::parse_argument( argc, argv, "--img", img_path );
//    float scale = 0.1f;
//    pcl::console::parse_argument( argc, argv, "--scale", scale );

//    return GF2::Solver::run( img_path, scale, {0, M_PI_2, M_PI}, argc, argv );
}

