#include "pearl/pearl.h"

#include  <vector>

#include "globfit2/globOpt_types.h"
#include "globfit2/util/parse.h"    // console::parse_argument
#include "globfit2/io/io.h"         // readPrimitives, readPoints
#include "globfit2/util/containers.hpp" // add
#include "globfit2/simple_types.h"


template int
am::Pearl::run(
        std::vector<int>                                        & labels
        , std::vector<GF2::LinePrimitive>      & lines
        , boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB> >            const& cloud
        , std::vector<int>   const* indices
        , am::Pearl::Params                    params
        , std::vector<std::vector<int> > *label_history
        , std::vector<std::vector<GF2::LinePrimitive> > *line_history
        );

template int
am::Pearl::run(
        std::vector<int>          & labels
        , std::vector<GF2::PlanePrimitive>      & lines
        , boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB> >            const& cloud
        , std::vector<int>   const* indices
        , am::Pearl::Params                    params
        , std::vector<std::vector<int> > *label_history
        , std::vector<std::vector<GF2::PlanePrimitive> > *line_history
        );

template < typename _PointContainerT
         , class    _InnerPrimitiveContainerT
         , typename _PrimitiveContainerT>
int pearlCli( int argc, char **argv )
{
    typedef typename _InnerPrimitiveContainerT::value_type    PrimitiveT;
    typedef typename _PointContainerT::value_type             PointPrimitiveT;
    typedef          std::map<GF2::GidT, _InnerPrimitiveContainerT> PrimitiveMapT;
    typedef          pcl::PointCloud<pcl::PointNormal>        PclCloudT;

    bool valid_input = true;

    std::string cloud_path       = "./cloud.ply",
                input_prims_path = "./patches.csv";

    am::Pearl::Params params;

    // parse
    {
        if (    (GF2::console::parse_argument( argc, argv, "--cloud", cloud_path) < 0)
             && (!boost::filesystem::exists(cloud_path)) )
        {
            std::cerr << "[" << __func__ << "]: " << "--cloud does not exist: " << cloud_path << std::endl;
            valid_input = false;
        }

        // primitives
        if (    (GF2::console::parse_argument( argc, argv, "-p"     , input_prims_path) < 0)
             && (GF2::console::parse_argument( argc, argv, "--prims", input_prims_path) < 0)
             && (!boost::filesystem::exists(input_prims_path)) )
        {
            std::cerr << "[" << __func__ << "]: " << "-p or --prims is compulsory" << std::endl;
            valid_input = false;
        }

        // scale
        if (    (GF2::console::parse_argument( argc, argv, "--scale", params.scale) < 0)
             && (GF2::console::parse_argument( argc, argv, "-sc"    , params.scale) < 0) )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale is compulsory" << std::endl;
            valid_input = false;
        }

        // weights
        pcl::console::parse_argument( argc, argv, "--unary", params.lambdas(0) );
        pcl::console::parse_argument( argc, argv, "--cmp"  , params.beta );
        valid_input &= pcl::console::parse_argument( argc, argv, "--pw"  , params.lambdas(2) ) >= 0;

        if (     !valid_input
              || (GF2::console::find_switch(argc,argv,"-h"    ))
              || (GF2::console::find_switch(argc,argv,"--help")) )
        {
            std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << "\n"
                      << "\t --cloud " << cloud_path << "\n"
                      << "\t -p,--prims " << input_prims_path << "\n"
                      << "\t -sc,--scale " << params.scale << "\n"
                      << "\t --pw " << params.lambdas(2) << "\n"
                      << "\t --cmp " << params.beta << "\n"
                      << "\n\t Example: ../pearl --scale 0.03 --cloud cloud.ply -p patches.csv --pw 1000 -cmp 1000\n"
                      << "\n";

            return EXIT_FAILURE;
        }
    } //...parse

    int err = EXIT_SUCCESS;

    // READ
    _PointContainerT points;
    PclCloudT::Ptr   pcl_cloud( new PclCloudT );
    if ( EXIT_SUCCESS == err )
    {
        err = GF2::io::readPoints<PointPrimitiveT>( points, cloud_path, &pcl_cloud );
        if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
    } //...read points

    _PrimitiveContainerT initial_primitives;
    PrimitiveMapT        patches;
    if ( EXIT_SUCCESS == err )
    {
        std::cout << "[" << __func__ << "]: " << "reading primitives from " << input_prims_path << "...";
        GF2::io::readPrimitives<PrimitiveT, _InnerPrimitiveContainerT>( initial_primitives, input_prims_path, &patches );
        std::cout << "reading primitives ok (#: " << initial_primitives.size() << ")\n";
    } //...read primitives

    // WORK
    {
        std::vector<int>          labels;
        _InnerPrimitiveContainerT primitives;
        err = am::Pearl::run( /* [out]        labels: */ labels
                            , /* [out]    primitives: */ primitives
                            , /* [in]      pcl_cloud: */ pcl_cloud
                            , /* [in]      p_indices: */ NULL
                            , /* [in]    pearlParams: */ params
                            , /* [out] label_history: */ NULL
                            , /* [out]  prim_history: */ (std::vector<std::vector<PrimitiveT> >*) NULL
                            //, /* [in]        patches: */ &initial_primitives
                            );

        PrimitiveMapT out_prims;
        std::cout<<"labels:";for(size_t vi=0;vi!=labels.size();++vi)std::cout<<labels[vi]<<" ";std::cout << "\n";
        for ( int pid = 0; pid != primitives.size(); ++pid )
        {
            const int gid = labels[pid];
            // assign point
            points[pid].setTag( PointPrimitiveT::TAGS::GID, gid );

            // add primitive
            if ( out_prims.find(gid) == out_prims.end() )
            {
                std::cout << "[" << __func__ << "]: " << "adding primitive[" << gid << "]: " << primitives[gid].toString() << std::endl;
                GF2::containers::add( out_prims, gid, primitives[gid] )
                        .setTag( PrimitiveT::TAGS::GID    , gid )
                        .setTag( PrimitiveT::TAGS::DIR_GID, gid );

            }
        }

        GF2::io::writeAssociations<PointPrimitiveT>( points, "./points_primitives.pearl.csv" );
        GF2::io::savePrimitives<PrimitiveT,typename _InnerPrimitiveContainerT::const_iterator>( out_prims, "./primitives.pearl.csv" );
    } //...work

    return err;
}


int main(int argc, char *argv[])
{
    std::cout << "hello pearl\n";
    if ( GF2::console::find_switch(argc,argv,"--3D") )
    {
        std::cout << "running planes" << std::endl;
        return pearlCli< GF2::PointContainerT
                       , GF2::_3d::InnerPrimitiveContainerT
                       , GF2::_3d::PrimitiveContainerT
                       >( argc, argv );
        return EXIT_FAILURE;
    }
    else
    {
        return pearlCli< GF2::PointContainerT
                       , GF2::_2d::InnerPrimitiveContainerT
                       , GF2::_2d::PrimitiveContainerT
                       >( argc, argv );
    } // if 3D

    return EXIT_FAILURE;
}
