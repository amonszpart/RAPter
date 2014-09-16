#include "globfit2/util/parse.h"

#include "globfit2/primitives/linePrimitive2.h"

namespace correspondance
{
    template <typename _PrimitiveT>
    int correspCli( int argc, char**argv )
    {
        // print usage
        if ( GF2::console::find_switch(argc,argv,"-h") || GF2::console::find_switch(argc,argv,"--help") )
        {
            std::cout << "Usage: "
                      << argv[0] << " --corresp"
                      << " --gt gt_prims.csv \n"
                      << " --p prims.csv"
                      << std::endl;
        } //...print usage

        // parse input
        std::string gt_path, prims_path;
        {
            if ( GF2::console::parse_argument(argc,argv,"--gt",gt_path) < 0 )
            {
                std::cerr << "[" << __func__ << "]: " << "need --gt gt_prims.csv to work" << std::endl;
                return EXIT_FAILURE;
            }

            if ( GF2::console::parse_argument(argc,argv,"-p",prims_path) < 0 )
            {
                std::cerr << "[" << __func__ << "]: " << "need -p gt_prims.csv to work" << std::endl;
                return EXIT_FAILURE;
            }
        } //...parse input

        {

        }

        return EXIT_SUCCESS;
    } //...correspCli()
} //...namespace correspondance

int corresp( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--corresp") )
    {
        return correspondance::correspCli<GF2::LinePrimitive2>(argc,argv);
    }

    return EXIT_FAILURE;
} //...corresp()
