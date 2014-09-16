#include "globfit2/util/parse.h"
#include "globfit2/globOpt_types.h" // _2d::, _3d::

// ________________correspondance___________________

#include "globfit2/io/io.h"         // readPrimitives()

namespace GF2 {
namespace correspondence
{
    template <typename Scalar, class _PrimitiveT, class _PointContainerT>
    inline Scalar estimateDistance( _PrimitiveT const& prim, _PrimitiveT const& gt_prim, _PointContainerT const& points, int pnt_gid, int gt_pnt_gid )
    {
        std::cout << "pos: " << prim.template pos().transpose() << ", gtpos: " << gt_prim.template pos().transpose() << ", norm: " << (prim.template pos() - gt_prim.template pos()).norm() << std::endl;
        return (prim.template pos() - gt_prim.template pos()).norm();
    }

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

        typedef std::pair< int   , int    > GidLid;   //!< Uniquely identifies an entry in PrimitiveMapT.
        typedef std::map < GidLid, GidLid > CorrespT; //!< Output type, contains primitive-gt correspondances. Key: <PrimitiveGid,PrimitiveLid> => Value: <GTGid,GTLid>

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
            _PrimitiveContainerT gt_primitives, primitives; // unused, so local scope
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

        // WORK
        {
            typedef typename PrimitiveMapT::const_iterator             outer_const_iterator;
            typedef typename _InnerPrimitiveContainerT::const_iterator inner_const_iterator;

            typedef typename _PointPrimitiveT::Scalar Scalar;
            typedef std::pair< GidLid , GidLid >      CostKey; // < primGid          ,primLid       > => <gtGid,gtLid>
            typedef std::map < CostKey, Scalar >      CostMap; // < <primGid,primLid>,<gtGid,gtLid> > => cost           // watch the order! <Prim, GT>
            CostMap costs;

            // calc distances
            int gid, gid_gt, lid, lid_gt;
            for ( outer_const_iterator outer_it0 = prims_map.begin(); outer_it0 != prims_map.end(); ++outer_it0 )
            {
                gid = (*outer_it0).first;
                lid = 0;
                for ( inner_const_iterator inner_it0 = (*outer_it0).second.begin(); inner_it0 != (*outer_it0).second.end(); ++inner_it0, ++lid )
                {
                    for ( outer_const_iterator outer_it1 = gt_prims_map.begin(); outer_it1 != gt_prims_map.end(); ++outer_it1 )
                    {
                        gid_gt = (*outer_it1).first;
                        lid_gt = 0;
                        for ( inner_const_iterator inner_it1 = (*outer_it1).second.begin(); inner_it1 != (*outer_it1).second.end(); ++inner_it1, ++lid_gt )
                        {
                            std::cout << "checking " << gid << "." << lid << " vs " << gid_gt << "." << lid_gt;
                            costs[ CostKey(GidLid(gid,lid),GidLid(gid_gt,lid_gt)) ] = estimateDistance<float>( *inner_it0, *inner_it1, points, PNT_GID, GT_PNT_GID );
                            std::cout << ": " << costs[ CostKey(GidLid(gid,lid),GidLid(gid_gt,lid_gt)) ] << std::endl;
                        } //...for gt_prims
                    } //...for gt_prims_map
                } //...for prims
            } //...for prims_map

            CorrespT corresps; // map<gidlid,gidlid>
            {
                // copy to list
                typedef std::pair<Scalar,CostKey> CostEntry;
                std::vector< CostEntry > cost_list;
                for ( typename CostMap::const_iterator it = costs.begin(); it != costs.end(); ++it )
                {
                    cost_list.push_back( CostEntry(it->second,it->first) );
                }

                // sort
                std::sort( cost_list.begin(), cost_list.end() );

                // select
                std::set<GidLid> taken_prims, taken_gts;
                for ( int i = 0; i < cost_list.size(); ++i )
                {
                    CostKey const& costKey    = cost_list[i].second;
                    GidLid  const& primGidLid = costKey.first;
                    GidLid  const& gtGidLid   = costKey.second;
                    std::cout << "cost_list[" << i << "]: " << cost_list[i].first << " for "
                              << primGidLid.first << "." << primGidLid.second << " - "
                              << gtGidLid.first << "." << gtGidLid.second << std::endl;

                    if ( (taken_prims.find(primGidLid) == taken_prims.end()) && (taken_gts.find(gtGidLid) == taken_gts.end()) )
                    {
                        std::cout << "chose " << cost_list[i].first << " for "
                                  << primGidLid.first << "." << primGidLid.second << " - "
                                  << gtGidLid.first << "." << gtGidLid.second << std::endl;

                        taken_prims.insert( primGidLid );
                        taken_gts  .insert( gtGidLid   );

                        if ( corresps.find(primGidLid) != corresps.end() )
                            std::cerr << "[" << __func__ << "]: " << "duplicate choice...should not happen!" << std::endl;

                        corresps[ primGidLid ] = gtGidLid;
                    }
                }
            }

            // print
            {
                for ( CorrespT::const_iterator it = corresps.begin(); it != corresps.end(); ++it )
                {
                    GidLid const& primGidLid = (*it).first;
                    GidLid const& gtGidLid = (*it).second;
                    std::cout << "prims[" << primGidLid.first << "][" << primGidLid.second << "]: " << prims_map[primGidLid.first][primGidLid.second].toString()
                              << " <--> "
                              << "gt[" << gtGidLid.first << "][" << gtGidLid.second << "]: " << gt_prims_map[gtGidLid.first][gtGidLid.second].toString()
                              << " with cost " << costs[ CostKey(primGidLid,gtGidLid) ]
                              << std::endl;
                }
            }
        } //...work

        return EXIT_SUCCESS;
    } //...correspCli()

} //...namespace correspondance
} //...namespace GF2

int corresp( int argc, char** argv )
{
    if ( GF2::console::find_switch(argc,argv,"--corresp") )
    {
        return GF2::correspondence::correspCli< GF2::_2d::PrimitiveT
                                              , GF2::_2d::InnerPrimitiveContainerT
                                              , GF2::_2d::PrimitiveContainerT
                                              , GF2::PointPrimitiveT
                                              , GF2::PointContainerT
                                              >( argc, argv );
    }

    return EXIT_FAILURE;
} //...corresp()
