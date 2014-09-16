#include "globfit2/util/parse.h"
#include "globfit2/globOpt_types.h" // _2d::, _3d::

// ________________correspondance___________________

#include "globfit2/io/io.h"         // readPrimitives()

namespace GF2 {
namespace correspondance
{
    template < typename _PrimitiveT
             , class    _InnerPrimitiveContainerT
             , class    _PrimitiveContainerT
             , class    _PointPrimitiveT
             , class    _PointContainerT
             >
    int correspCli( int argc, char**argv )
    {
        typedef std::map<int, _InnerPrimitiveContainerT> PrimitiveMapT;
        enum { PNT_GID    = _PointPrimitiveT::GID
             , GT_PNT_GID = _PointPrimitiveT::USER_ID1
             };

        int err = EXIT_SUCCESS;

        // print usage
        if ( GF2::console::find_switch(argc,argv,"-h") || GF2::console::find_switch(argc,argv,"--help") )
        {
            std::cout << "Usage: "
                      << argv[0] << " --corresp"
                      << " --gt gt_prims.csv \n"
                      << " --gta gt_points_primitives.csv\n"
                      << " --p prims.csv\n"
                      << " --pa points_primitives.csv\n"
                      << " --cloud cloud.ply\n"
                      ;
            return err;
        } //...print usage

        // parse input
        std::string gt_path,
                    prims_path,
                    cloud_path = "./cloud.ply",
                    assoc_path,
                    gt_assoc_path;
        {
            if ( GF2::console::parse_argument(argc,argv,"--gt",gt_path) < 0 )
            {
                std::cerr << "[" << __func__ << "]: " << "need --gt gt_prims.csv to work" << std::endl;
                return EXIT_FAILURE;
            }

            if ( GF2::console::parse_argument(argc,argv,"--p",prims_path) < 0 )
            {
                std::cerr << "[" << __func__ << "]: " << "need --p gt_prims.csv to work" << std::endl;
                return EXIT_FAILURE;
            }

            if ( GF2::console::parse_argument(argc,argv,"--gta",gt_assoc_path) < 0 )
            {
                std::cerr << "[" << __func__ << "]: " << "need --gta gt_points_primitives.csv to work" << std::endl;
                return EXIT_FAILURE;
            }

            if ( GF2::console::parse_argument(argc,argv,"--pa", assoc_path) < 0 )
            {
                std::cerr << "[" << __func__ << "]: " << "need --pa points_primitives.csv to work" << std::endl;
                return EXIT_FAILURE;
            }

            if ( (GF2::console::parse_argument(argc,argv,"--cloud",cloud_path) < 0) && !boost::filesystem::exists(cloud_path) )
            {
                std::cerr << "[" << __func__ << "]: " << "need --cloud cloud.ply to work" << std::endl;
                return EXIT_FAILURE;
            }
        } //...parse input

        _PointContainerT     points; // need two versions to hold the associations...:-S
        _PrimitiveContainerT gt_primitives, primitives;
        PrimitiveMapT        gt_prims_map, prims_map;
        // read input
        {
            // Read points
            if ( EXIT_SUCCESS == err )
            {
                err = io::readPoints<_PointPrimitiveT>( points, cloud_path );
                if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
            } //...read points

            std::vector<std::pair<int,int> > points_primitives;

            // read GT associations
            {
                io::readAssociations( points_primitives, gt_assoc_path, NULL );
                for ( size_t i = 0; i != points.size(); ++i )
                {
                    // store association in point
                    points[i].setTag( GT_PNT_GID, points_primitives[i].first );
                }
            }

            // read associations
            {
                io::readAssociations( points_primitives, assoc_path, NULL );
                for ( size_t i = 0; i != points.size(); ++i )
                {
                    // store association in point
                    points[i].setTag( PNT_GID, points_primitives[i].first );
                }
            }

            // read primitives
            {
                // GT
                std::cout << "[" << __func__ << "]: " << "reading GT primitives from " << gt_path << "...";
                io::readPrimitives<_PrimitiveT, _InnerPrimitiveContainerT>( gt_primitives, gt_path, &gt_prims_map );
                std::cout << "reading GT primitives ok (#: " << gt_prims_map.size() << ")\n";

                // solution
                std::cout << "[" << __func__ << "]: " << "reading primitives from " << prims_path << "...";
                io::readPrimitives<_PrimitiveT, _InnerPrimitiveContainerT>( primitives, prims_path, &prims_map );
                std::cout << "reading primitives ok (#: " << prims_map.size() << ")\n";
            } //...read primitives
        }

        return EXIT_SUCCESS;
    } //...correspCli()

} //...namespace correspondance
} //...namespace GF2

int corresp( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--corresp") )
    {
        return GF2::correspondance::correspCli< GF2::_2d::PrimitiveT
                                              , GF2::_2d::InnerPrimitiveContainerT
                                              , GF2::_2d::PrimitiveContainerT
                                              , GF2::PointPrimitiveT
                                              , GF2::PointContainerT
                                              >( argc, argv );
    }

    return EXIT_FAILURE;
} //...corresp()
