#include <iostream>

#include <pcl/console/parse.h>

#include "optimization/qp/solver.h"
#include "globfit2/optimization/problemSetup.h"

#include "globfit2/io/io.h"

int segment3D( int argc, char** argv );
int segment2D( int argc, char** argv );
int subsample( int argc, char** argv );
int merge    ( int argc, char** argv );

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
                  << "\t--merge3D\n"
                  << "\t--gfit\n"
                  << "\t--show\n"
                  << std::endl;

        return EXIT_SUCCESS;
    }
#if GF2_WITH_SAMPLE_INPUT
    else if ( pcl::console::find_switch(argc,argv,"--sample-input") )
    {
        return GF2::Solver::sampleInput( argc, argv );
    }
#endif
    else if ( pcl::console::find_switch(argc,argv,"--segment") )
    {
       return segment2D( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--segment3D") )
    {
        return segment3D( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--generate") )
    {
        return GF2::Solver::generateCli( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--formulate") )
    {
        return GF2::ProblemSetup::formulateCli<GF2::Solver::PrimitiveContainerT, GF2::Solver::PointContainerT>( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--solver") )
    {
        return GF2::Solver::solve( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--gfit") )
    {
        return GF2::Solver::datafit( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--merge") || pcl::console::find_switch(argc,argv,"--merge3D") )
    {
        return merge(argc, argv);
    }
    else if ( pcl::console::find_switch(argc,argv,"--show") )
    {
        std::cerr << "[" << __func__ << "]: " << "the show option has been moved to a separate executable, please use thatt one" << std::endl;
        return 1;
        //return GF2::Solver::show( argc, argv );
    }
    else if ( pcl::console::find_switch(argc,argv,"--subsample") )
    {
        return subsample( argc, argv );
    }

    return 1;

    // --show --dir . --cloud cloud.ply --scale 0.05f --assoc points_primitives.txt --use-tags --no-clusters --prims primitives.bonmin.txt

//    std::string img_path( "input2.png" );
//    pcl::console::parse_argument( argc, argv, "--img", img_path );
//    float scale = 0.1f;
//    pcl::console::parse_argument( argc, argv, "--scale", scale );

//    return GF2::Solver::run( img_path, scale, {0, M_PI_2, M_PI}, argc, argv );
}

