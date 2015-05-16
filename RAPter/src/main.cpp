#include <iostream>

#include "rapter/util/parse.h"

int subsample ( int argc, char** argv ); // subsample.cpp
int segment   ( int argc, char** argv ); // segment.cpp
int generate  ( int argc, char** argv ); // generate.cpp
int generate3D( int argc, char** argv ); // generate3D.cpp
int formulate ( int argc, char** argv ); // problemSetup.cpp
int formulate3D( int argc, char** argv ); // problemSetup.cpp
int solve     ( int argc, char** argv ); // solve.cpp
int solve3D   ( int argc, char** argv ); // solve3D.cpp
int merge     ( int argc, char** argv ); // merge.cpp
//int datafit   ( int argc, char** argv ); // datafit.cpp
//int reassign  ( int argc, char** argv );
int represent ( int argc, char** argv ); // represent.cpp

int main( int argc, char *argv[] )
{
    if ( (argc == 2) &&
         (   (rapter::console::find_switch(argc,argv,"--help"))
          || (rapter::console::find_switch(argc,argv,"-h"    ))
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
                  << "\t--corresp\n"
                  << "\t--represent[3D]"
                  //<< "\t--show\n"
                  << std::endl;

        return EXIT_SUCCESS;
    }
    else if ( rapter::console::find_switch(argc,argv,"--segment") || rapter::console::find_switch(argc,argv,"--segment3D") )
    {
       return segment( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--generate") )
    {
        return generate(argc,argv);
    }
    else if ( rapter::console::find_switch(argc,argv,"--generate3D") )
    {
        return generate3D(argc,argv);
    }
    else if ( rapter::console::find_switch(argc,argv,"--formulate") )
    {
        return formulate( argc, argv );
        //return rapter::ProblemSetup::formulateCli<rapter::Solver::PrimitiveContainerT, rapter::Solver::PointContainerT>( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--formulate3D") )
    {
        return formulate3D( argc, argv );
        //return rapter::ProblemSetup::formulateCli<rapter::Solver::PrimitiveContainerT, rapter::Solver::PointContainerT>( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--solver") ) // Note: "solver", not "solve" :-S
    {
        return solve( argc, argv );
        //return rapter::Solver::solve( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--solver3D") ) // Note: "solver", not "solve" :-S
    {
        return solve3D( argc, argv );
    }
//    else if ( rapter::console::find_switch(argc,argv,"--datafit") || rapter::console::find_switch(argc,argv,"--datafit3D") )
//    {
//        return datafit( argc, argv );
//        //return rapter::Solver::datafit( argc, argv );
//    }
    else if ( rapter::console::find_switch(argc,argv,"--merge") || rapter::console::find_switch(argc,argv,"--merge3D") )
    {
        return merge(argc, argv);
    }
    else if ( rapter::console::find_switch(argc,argv,"--show") )
    {
        std::cerr << "[" << __func__ << "]: " << "the show option has been moved to a separate executable, please use thatt one" << std::endl;
        return 1;
        //return rapter::Solver::show( argc, argv );
    }
    else if ( rapter::console::find_switch(argc,argv,"--subsample") )
    {
        return subsample( argc, argv );
    }
//    else if ( rapter::console::find_switch(argc,argv,"--reassign") )
//    {
//        return reassign( argc, argv );
//    }
    else if ( rapter::console::find_switch(argc,argv,"--represent") || rapter::console::find_switch(argc,argv,"--represent3D")
              || rapter::console::find_switch(argc,argv,"--representBack") || rapter::console::find_switch(argc,argv,"--representBack3D") )
    {
        return represent( argc, argv );
    }
//    else if ( rapter::console::find_switch(argc,argv,"--corresp") || rapter::console::find_switch(argc,argv,"--corresp3D") )
//    {
//        return corresp( argc, argv );
//    }

    std::cerr << "[" << __func__ << "]: " << "unrecognized option" << std::endl;
    return 1;

    // --show --dir . --cloud cloud.ply --scale 0.05f --assoc points_primitives.txt --use-tags --no-clusters --prims primitives.bonmin.txt

//    std::string img_path( "input2.png" );
//    pcl::console::parse_argument( argc, argv, "--img", img_path );
//    float scale = 0.1f;
//    pcl::console::parse_argument( argc, argv, "--scale", scale );

//    return rapter::Solver::run( img_path, scale, {0, M_PI_2, M_PI}, argc, argv );
}

