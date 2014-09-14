#include <iostream>

//#include <pcl/console/parse.h>
#include "globfit2/util/parse.h"

//#include "optimization/qp/solver.h"
//#include "globfit2/optimization/problemSetup.h"


int segment  ( int argc, char** argv ); // segment3D.cpp //TODO: rename
int subsample( int argc, char** argv ); // subsample.cpp
int merge    ( int argc, char** argv ); // merge.cpp
int generate ( int argc, char** argv ); // generate.cpp
int formulate( int argc, char** argv ); // problemSetup.cpp
int solve    ( int argc, char** argv ); // solve.cpp
int datafit  ( int argc, char** argv ); // datafit.cpp

int main( int argc, char *argv[] )
{
    if ( (argc == 2) &&
         (   (GF2::console::find_switch(argc,argv,"--help"))
          || (GF2::console::find_switch(argc,argv,"-h"    ))
         )
       )
    {
        std::cout << "[Usage]:\n"
                  << "\t--generate\n"
                  << "\t--generate3D\n"
                  << "\t--formulate\n"
                  << "\t--formulate3D\n"
                  << "\t--solver mosek|bonmin|gurobi\n"
                  << "\t--solver3D mosek|bonmin|gurobi\n"
                  << "\t--merge\n"
                  << "\t--merge3D\n"
                  << "\t--datafit\n"
                  //<< "\t--show\n"
                  << std::endl;

        return EXIT_SUCCESS;
    }
    else if ( GF2::console::find_switch(argc,argv,"--segment") || GF2::console::find_switch(argc,argv,"--segment3D") )
    {
       return segment( argc, argv );
    }
    else if ( GF2::console::find_switch(argc,argv,"--generate") || GF2::console::find_switch(argc,argv,"--generate3D") )
    {
        return generate(argc,argv);
    }
    else if ( GF2::console::find_switch(argc,argv,"--formulate") || GF2::console::find_switch(argc,argv,"--formulate3D"))
    {
        return formulate( argc, argv );
        //return GF2::ProblemSetup::formulateCli<GF2::Solver::PrimitiveContainerT, GF2::Solver::PointContainerT>( argc, argv );
    }
    else if ( GF2::console::find_switch(argc,argv,"--solver") || GF2::console::find_switch(argc,argv,"--solver3D") ) // Note: "solver", not "solve" :-S
    {
        return solve( argc, argv );
        //return GF2::Solver::solve( argc, argv );
    }
    else if ( GF2::console::find_switch(argc,argv,"--datafit") || GF2::console::find_switch(argc,argv,"--datafit3D") )
    {
        return datafit( argc, argv );
        //return GF2::Solver::datafit( argc, argv );
    }
    else if ( GF2::console::find_switch(argc,argv,"--merge") || GF2::console::find_switch(argc,argv,"--merge3D") )
    {
        return merge(argc, argv);
    }
    else if ( GF2::console::find_switch(argc,argv,"--show") )
    {
        std::cerr << "[" << __func__ << "]: " << "the show option has been moved to a separate executable, please use thatt one" << std::endl;
        return 1;
        //return GF2::Solver::show( argc, argv );
    }
    else if ( GF2::console::find_switch(argc,argv,"--subsample") )
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

