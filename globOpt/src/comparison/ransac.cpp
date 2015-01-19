#include  <vector>

#include "globfit2/globOpt_types.h"
#include "globfit2/util/parse.h"    // console::parse_argument
#include "globfit2/io/io.h"         // readPrimitives, readPoints
#include "globfit2/util/containers.hpp" // add
#include "globfit2/processing/util.hpp" // getpop
#include "schnabelEnv.h"
#include "../../src/schnabelEnv.cpp"

#if GF2_USE_PCL
    #include "pcl/sample_consensus/model_types.h" // sacmodel
    #include "pcl/sample_consensus/sac_model_line.h"
    #include "pcl/sample_consensus/sac_model_plane.h"
    #include "pcl/sample_consensus/ransac.h"
    #include "pcl/sample_consensus/sac_model.h"
#endif // GF2_USE_PCL

template <typename _Scalar>
struct RansacParams
{
    enum ALGO { RANSAC };

    ALGO            algo = RANSAC;
    _Scalar         scale;
    pcl::SacModel   modelType;
};

#include "globfit2/io/inputParser.hpp"
#include "globfit2/comparison/impl/reassign.hpp"

template < typename _PointContainerT
         , class    _InnerPrimitiveContainerT
         , typename _PrimitiveContainerT
         , class    _PclCloudT >
int ransacCli( int argc, char **argv )
{
    typedef typename _InnerPrimitiveContainerT::value_type    PrimitiveT;
    typedef typename _PointContainerT::value_type             PointPrimitiveT;
    typedef          std::map<GF2::GidT, _InnerPrimitiveContainerT> PrimitiveMapT;
    typedef typename PointPrimitiveT::Scalar                  Scalar;
    typedef typename _PclCloudT::PointType                    PclPointT;

    bool valid_input = true;
    bool usePatches = false;
    int maxPlanes = 0;

//    std::string cloud_path       = "./cloud.ply",
//                input_prims_path = "./patches.csv",
//                associations_path = "./points_primitives.csv";

    RansacParams<Scalar>    params;
    std::string             algo_str = "ransac";

    _PointContainerT         points;
    typename _PclCloudT::Ptr pclCloud;
    _PrimitiveContainerT     initial_primitives;
    PrimitiveMapT            patches;

    // parse
    {
        valid_input = EXIT_SUCCESS == GF2::parseInput<_InnerPrimitiveContainerT,_PclCloudT>( points, pclCloud, initial_primitives, patches, params, argc, argv );
        params.modelType = pcl::console::find_switch( argc, argv, "--3D" ) ? pcl::SacModel::SACMODEL_PLANE
                                                                           : pcl::SacModel::SACMODEL_LINE;
        usePatches = pcl::console::find_switch( argc, argv, "--use-patches");

        // weights
        //pcl::console::parse_argument( argc, argv, "--unary", params.lambdas(0) );
        //pcl::console::parse_argument( argc, argv, "--cmp"  , params.beta );
        //valid_input &= pcl::console::parse_argument( argc, argv, "--pw"  , params.lambdas(2) ) >= 0;

        pcl::console::parse_argument( argc, argv, "--algo", algo_str );
        if ( !algo_str.compare("ransac") )
            params.algo = RansacParams<Scalar>::RANSAC;
        else
            std::cerr << "[" << __func__ << "]: " << "unrecognized ransac algorithm" << std::endl;

        pcl::console::parse_argument( argc, argv, "--max-prims", maxPlanes );

        if (     !valid_input
              || (GF2::console::find_switch(argc,argv,"-h"    ))
              || (GF2::console::find_switch(argc,argv,"--help")) )
        {
            std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << "\n"
                      << "\t --cloud " << /*cloud_path <<*/ "\n"
                      << "\t -p,--prims " << /*input_prims_path <<*/ "\n"
                      << "\t -a,--assoc " << /*associations_path <<*/ "\n"
                      << "\t -sc,--scale " << params.scale << "\n"
                      << "\t[ --use-patches " << (usePatches?"YES":"NO") << "]\t Use inital segmentation\n"
                      << "\t[--max-prims " << maxPlanes << "]\t Limit #output primitives"
                      << "\t --3D\n"
                      //<< "\t --pw " << params.lambdas(2) << "\n"
                      //<< "\t --cmp " << params.beta << "\n"
                      << "\t Example: ../ransac --scale 0.03 --cloud cloud.ply -p patches.csv -a points_primitives.csv"
                      << "\n";

            return EXIT_FAILURE;
        }
    } //...parse

    int err = EXIT_SUCCESS;

    PrimitiveMapT out_prims;
    std::vector<int> inliers;

    // created RandomSampleConsensus object and compute the appropriated model
    GF2::GidPidVectorMap populations;
    GF2::processing::getPopulations( populations, points );

    if ( !usePatches )
    {
        for ( size_t pid = 0; pid != points.size(); ++pid )
            points[pid].setTag( PointPrimitiveT::TAGS::GID, PointPrimitiveT::TAG_UNSET );
        if ( params.modelType == pcl::SacModel::SACMODEL_PLANE
             || params.modelType == pcl::SacModel::SACMODEL_LINE )
        {
            if ( params.algo == RansacParams<Scalar>::RANSAC )
            {
//                typename _PclCloudT::Ptr restCloud( new _PclCloudT() );
//                pcl::copyPointCloud( pclCloud, restCloud );

                //out_prims
                std::vector<int> indices( pclCloud->size() );
                for ( size_t pid = 0; pid != pclCloud->size(); ++pid )
                    indices[pid] = pid;
                while ( (indices.size() > 2) && ((maxPlanes == 0) || (out_prims.size() < maxPlanes)) )
                {
//                    found = false;
                    //typename pcl::SampleConsensusModelLine<PclPointT>::Ptr model( new pcl::SampleConsensusModelLine<PclPointT> (pclCloud, indices) );
                    //typename pcl::SampleConsensusModelPlane<PclPointT>::Ptr model( new pcl::SampleConsensusModelPlane<PclPointT> (pclCloud, indices) );
                    typename pcl::SampleConsensusModel<PclPointT>::Ptr model;
                    if ( params.modelType == pcl::SacModel::SACMODEL_PLANE )
                        model.reset( new pcl::SampleConsensusModelPlane<PclPointT> (pclCloud, indices) );
                    else // LINE
                        model.reset( new pcl::SampleConsensusModelLine<PclPointT> (pclCloud, indices) );
                    pcl::RandomSampleConsensus<PclPointT> ransac( model );
                    ransac.setDistanceThreshold (params.scale);
                    if ( ransac.computeModel() )
                    {
                        ransac.getInliers(inliers);
                        // save line
                        if ( inliers.size() )
                        {
                            std::cout << "ret: " << ransac.model_coefficients_.transpose() << std::endl;

                            PrimitiveT* primToAdd = NULL;
                            if ( params.modelType == pcl::SacModel::SACMODEL_PLANE )
                            {
                                Eigen::Matrix<Scalar,3,1> nrm = ransac.model_coefficients_.template head<3>();
                                nrm.normalize();
                                primToAdd = new PrimitiveT( Eigen::Matrix<Scalar,3,1>::Zero() - nrm * ransac.model_coefficients_(3), nrm );
                            }
                            else
                            { //LINE
                                primToAdd = new PrimitiveT( ransac.model_coefficients_ );
                            }

                            int gid = out_prims.size();

                            GF2::containers::add( out_prims, gid, *primToAdd )
                                    .setTag( PrimitiveT::TAGS::GID    , gid )
                                    .setTag( PrimitiveT::TAGS::DIR_GID, gid );
                            std::cout << "created " << out_prims[gid].back().toString() << ", distO: " << out_prims[gid].back().getDistance( Eigen::Matrix<Scalar,3,1>::Zero() )  << std::endl;
//                            found = true;

                            // add inliers
                            auto indicesIt = indices.begin();
                            for ( auto inlierIt = inliers.begin(); inlierIt != inliers.end(); ++inlierIt )
                            {
                                // assign point
                                points.at( *inlierIt ).setTag( PointPrimitiveT::TAGS::GID, gid );
                                // delete from cloud indices
                                while ( (*indicesIt < *inlierIt) && (indicesIt != indices.end()) )
                                    ++indicesIt;
                                if ( *indicesIt == *inlierIt )
                                    indices.erase( indicesIt );
                            }
                            inliers.clear();

                            if ( primToAdd ) delete primToAdd;
                        } //...if inliners
                    }
                    else
                        inliers.clear();
                }
            }
            else
                std::cerr << "unimplemented algorithm " << algo_str << std::endl;
        }
    }
    else
    {
        for ( typename PrimitiveMapT::const_iterator patch_it = patches.begin(); patch_it != patches.end(); ++patch_it )
        {
            //for ( _InnerPrimitiveContainerT::const_iterator prim_it = patch_it->begin(); prim_it != patch_it->end(); ++prim_it )
            const int gid = patch_it->first;
            if ( populations[gid].size() )
            {
                if ( params.modelType == pcl::SacModel::SACMODEL_PLANE )
                {
                    if ( params.algo == RansacParams<Scalar>::RANSAC )
                    {
                        std::vector<int> indices;
                        std::copy( populations[gid].begin(), populations[gid].end(), indices.begin() );

                        typename pcl::SampleConsensusModelPlane<PclPointT>::Ptr model( new pcl::SampleConsensusModelPlane<PclPointT> (pclCloud, indices) );
                        pcl::RandomSampleConsensus<PclPointT> ransac (model);
                        ransac.setDistanceThreshold (params.scale);
                        if ( ransac.computeModel() )
                        {
                            ransac.getInliers(inliers);
                            // save line
                            if ( inliers.size() )
                            {
                                std::cout << "ret: " << ransac.model_coefficients_.transpose() << std::endl;

                                Eigen::Matrix<Scalar,3,1> nrm = ransac.model_coefficients_.template head<3>();
                                nrm.normalize();
                                GF2::containers::add( out_prims, gid, PrimitiveT( Eigen::Matrix<Scalar,3,1>::Zero() - nrm * ransac.model_coefficients_(3), nrm) )
                                        .setTag( PrimitiveT::TAGS::GID    , gid )
                                        .setTag( PrimitiveT::TAGS::DIR_GID, gid );
                                std::cout << "created " << out_prims[gid].back().toString() << ", distO: " << out_prims[gid].back().getDistance( Eigen::Matrix<Scalar,3,1>::Zero() )  << std::endl;
                            }
                        }
                        else
                            inliers.clear();
                    }
                    else
                        std::cerr << "unimplemented algorithm " << algo_str << std::endl;
                }
                else if (params.modelType == pcl::SacModel::SACMODEL_LINE )
                {
                    if ( params.algo == RansacParams<Scalar>::RANSAC )
                    {
                        std::vector<int> indices;
                        std::copy( populations[gid].begin(), populations[gid].end(), indices.begin() );

                        std::cout << "[" << __func__ << "]: " << "doing ransac lines" << std::endl;
                        std::cout<<"populations["<<gid<<"]:";for(size_t vi=0;vi!=populations[gid].size();++vi)std::cout<<populations[gid][vi]<<" ";std::cout << "\n";
                        typename pcl::SampleConsensusModelLine<PclPointT>::Ptr model( new pcl::SampleConsensusModelLine<PclPointT> (pclCloud, indices) );

                        pcl::RandomSampleConsensus<PclPointT> ransac( model );
                        ransac.setDistanceThreshold( params.scale );
                        if ( ransac.computeModel() )
                        {
                            ransac.getInliers( inliers );

                            // save line
                            if ( inliers.size() )
                            {
                                GF2::containers::add( out_prims, gid, PrimitiveT( ransac.model_coefficients_.template head<3>(), ransac.model_coefficients_.template segment<3>(3)) )
                                        .setTag( PrimitiveT::TAGS::GID    , gid )
                                        .setTag( PrimitiveT::TAGS::DIR_GID, gid );
                                std::cout << "created " << out_prims[gid].back().toString() << std::endl;
                            }
                        }
                        else
                            inliers.clear();
                    }
                    else
                        std::cerr << "unimplemented algorithm " << algo_str << std::endl;
                }

    #if 0
                // copies all inliers of the model computed to another PointCloud
                std::cout << "copying " << inliers.size() << " points" << std::endl;
                PclCloudT::Ptr out_cloud( new PclCloudT() );
                pcl::copyPointCloud<PclPointT>(*pcl_cloud, inliers, *out_cloud);

                pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
                vptr->setBackgroundColor( .5, .6, .6 );
                vptr->addPointCloud<PclPointT>( out_cloud );
                vptr->spin();
    #endif
            } // ...if popsize
        } //...for each patch
        } //...if use patch

    GF2::io::savePrimitives<PrimitiveT,typename _InnerPrimitiveContainerT::const_iterator>( out_prims, "./primitives.ransac.csv" );
    GF2::io::writeAssociations<PointPrimitiveT>( points, "./points_primitives.ransac.csv" );

    return err;
}

#include <chrono>
#define TIC auto start = std::chrono::system_clock::now();
#define RETIC start = std::chrono::system_clock::now();
#define TOC(title,it) { std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start; \
                     std::cout << title << ": " << elapsed_seconds.count()/it << " s" << std::endl; }

template < class _PrimitiveContainerT
         , class _PointContainerT
         , class _Scalar
//         , class _PrimitiveMapT
         >
inline int reassign( _PointContainerT &points, _PrimitiveContainerT const& primitives, _Scalar const scale/*, _PrimitiveMapT const& patches*/ )
{
    typedef typename _PointContainerT::value_type PointPrimitiveT;

#if 1
    std::cout << "starting assignment" << std::endl; fflush( stdout );

    // distances between a point patch (point.gid), and a new primitive (primitive.gid)
    std::map< std::pair<int,int>, std::pair<_Scalar,int> > dists; // < <point.gid, primitives.gid>, <sumdist,|points|> >

    TIC
    #pragma omp for schedule(dynamic,8)
    for ( int pid = 0; pid < points.size(); ++pid )
    {
        float min_dist = FLT_MAX, tmp;
        int   min_gid = 0;
        for ( int gid = 0; gid != primitives.size(); ++gid )
        {
            tmp = std::abs( primitives[gid].getDistance(points[pid].template pos()) );
            if ( tmp < 0.f )
                std::cerr << "asdf: " << tmp << std::endl;
            std::pair<int,int> key( points[pid].getTag(PointPrimitiveT::TAGS::GID), gid );
            #pragma omp critical
            {
                dists[ key ].first += tmp;
                dists[ key ].second++;
            }

            if ( (tmp < min_dist) && (tmp < scale) )
            {
                min_dist = tmp;
                min_gid = gid;
            }
        }
        points[pid].setTag( PointPrimitiveT::TAGS::GID, min_gid );
    }
    TOC("reassign omp",1)
    std::cout << "finishing assignment" << std::endl;
    #endif

    return 0;
}

template < typename _PointContainerT
         , class    _InnerPrimitiveContainerT
         , typename _PrimitiveContainerT
         , class    _PclCloudT >
inline int schnabelCli( int argc, char** argv )
{
    typedef typename _PclCloudT::PointType                    PclPointT;
    typedef          std::map<GF2::GidT, _InnerPrimitiveContainerT> PrimitiveMapT;
    typedef typename _PointContainerT::value_type             PointPrimitiveT;
    typedef typename _InnerPrimitiveContainerT::value_type   PrimitiveT;
    typedef typename _PclCloudT::PointType                    PclPointT;

    bool extrude2D = false;
    if ( GF2::console::find_switch(argc,argv,"--schnabel2D") )
        extrude2D = true;
    if ( extrude2D )
        std::cout << "hello Schnabel2D\n";
    else
        std::cout << "hello Schnabel3D\n";

    RansacParams<Scalar>     params;
    _PointContainerT         points, outPoints;
    typename _PclCloudT::Ptr pcl_cloud;
    _PrimitiveContainerT     initial_primitives;
    PrimitiveMapT            patches;
    int                      min_support_arg = 300;
    int pointMultiplier = 50;
    GF2::console::parse_argument( argc, argv, "--point-mult", pointMultiplier);

    // parse
    {
        bool valid_input = !GF2::parseInput<_InnerPrimitiveContainerT,_PclCloudT>( points, pcl_cloud, initial_primitives, patches, params, argc, argv );

        if ( (GF2::console::parse_argument( argc, argv, "--minsup", min_support_arg) < 0) )
        {
            std::cerr << "[" << __func__ << "]: " << "--minsup required (300?)" << std::endl;
            valid_input = false;
        }

        if (     !valid_input
              || (GF2::console::find_switch(argc,argv,"-h"    ))
              || (GF2::console::find_switch(argc,argv,"--help")) )
        {
            std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << "\n"
                      << "\t --cloud " << /*cloud_path <<*/ "\n"
                      << "\t -p,--prims " << /*input_prims_path <<*/ "\n"
                      << "\t -a,--assoc " << /*associations_path <<*/ "\n"
                      << "\t -sc,--scale " << params.scale << "\n"
                      << "\t --minsup " << min_support_arg << "\n"
                      << "\t --point-mult " << pointMultiplier << "]\t extrude2D add this many points for each input point\n"
                      << "\t Example: ../ransac --schnabel3D --scale 0.03 --cloud cloud.ply -p patches.csv -a points_primitives.csv"
                      << "\n";

            return EXIT_FAILURE;
        }
    }

    int err = EXIT_SUCCESS;

    // debug visualize
    pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
    {
        pcl::PointCloud<pcl::PointXYZ> tmp;
        for ( size_t pid = 0; pid != pcl_cloud->size(); ++pid )
        {
            tmp.push_back( pcl::PointXYZ( pcl_cloud->at(pid).x
                                          , pcl_cloud->at(pid).y
                                          , pcl_cloud->at(pid).z) );
        }
        vptr->setBackgroundColor( .5, .6, .6 );
        vptr->addPointCloud<pcl::PointXYZ>( tmp.makeShared(), "asdf");
        vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4.f );

        for ( size_t i = 0; i != std::min(50UL,points.size()); ++i )
        {
            std::cout << "\tpoints[" << i << "]: " << points[i].toString() << "\n\tcloud[" << i << "]: " << pcl_cloud->at(i).getVector3fMap().transpose() << std::endl;
        }

        vptr->spinOnce(100);
    }

    typedef std::map<GF2::PidT,GF2::GidT> PidGidT;
     std::vector<GF2::PlanePrimitive> planes;
     PidGidT                          pidGid;

     err = GF2::SchnabelEnv::run<pcl::PointCloud<PclPointT> >( planes
                               //, pidGid
                               , outPoints
                               , points
                               , pcl_cloud
                               , params.scale
                               , min_support_arg
                               , false
                               , extrude2D
                               , pointMultiplier
                               );
    PrimitiveMapT out_prims;
    for ( size_t gid = 0; gid != planes.size(); ++gid )
    {
        GF2::containers::add( out_prims, gid, planes[gid] )
                .setTag( PrimitiveT::TAGS::GID    , gid )
                .setTag( PrimitiveT::TAGS::DIR_GID, gid );
        std::cout << "added " << out_prims[gid].back().toString() << std::endl;

        if ( extrude2D )
        {

        }
    }



//    for ( PidGidT::const_iterator it = pidGid.begin(); it != pidGid.end(); ++it )
//    {
//        points.at( it->first ).setTag( PointPrimitiveT::TAGS::GID, it->second );
//    }
#if 0
    reassign( points, planes, params.scale );

    // debug assignments
    pcl::PointXYZ pnt, plane_pnt;
    for ( int pid = 0; pid != static_cast<int>(points.size()); ++pid )
    {
        if ( !(pid % 1000) )
        {
            char name[255];
            sprintf( name, "pointarrow%06d", pid );
            pnt.getVector3fMap() = points[pid].pos();
            const int gid = points[pid].getTag( PointPrimitiveT::TAGS::GID );
            plane_pnt.getVector3fMap() = out_prims[gid].at(0).pos();
            vptr->addArrow( plane_pnt, pnt,  .1, .9, .5, false, name );
        }
    }

//        for ( PidGidT::const_iterator it = pidGid.begin(); it != pidGid.end(); ++it )
//        {
//            points[it->first].setTag( PointPrimitiveT::TAGS::GID, it->second );
//        }
    for ( size_t i = 0; i != std::min(50UL,planes.size()); ++i )
    {
        char name[255];
        sprintf( name, "plane%lu", i );
        vptr->addPlane( *(planes[i].modelCoefficients()), planes[i]. pos()(0), planes[i]. pos()(1), planes[i]. pos()(2), name );
    }
    vptr->addCoordinateSystem(0.5);
    vptr->spin();
#endif

    char path[ 256 ];
    sprintf( path, "./schnabel_minsup%d.primitives.csv", min_support_arg );
    std::cout << " writing " << out_prims.size() << " primitives to " << path << std::endl;
    GF2::io::savePrimitives<PrimitiveT,std::vector<GF2::PlanePrimitive>::const_iterator >( out_prims, std::string(path) );

    sprintf( path, "./schnabel_minsup%d.points_primitives.csv", min_support_arg );
    std::cout << " writing " << outPoints.size() << " assignments to " << path << std::endl;
    GF2::io::writeAssociations<PointPrimitiveT>( outPoints, path );

    sprintf( path, "./schnabel_minsup%d.cloud.ply", min_support_arg );
    std::cout << " writing " << outPoints.size() << " points to " << path << std::endl;
    GF2::io::writePoints<PointPrimitiveT>( outPoints, path );

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
    typedef          pcl::PointNormal                         PclPointT;
    typedef          pcl::PointCloud<PclPointT>               PclCloudT;

    std::cout << "hello ransac\n";
    if ( GF2::console::find_switch(argc,argv,"--schnabel3D") || GF2::console::find_switch(argc,argv,"--schnabel2D") )
    {
        return schnabelCli< GF2::PointContainerT
                          , GF2::_3d::InnerPrimitiveContainerT
                          , GF2::_3d::PrimitiveContainerT
                          , PclCloudT
                          >( argc, argv);
    }
    else if ( GF2::console::find_switch(argc,argv,"--assign") )
    {
        return GF2::reassignCli< GF2::PointContainerT
                               , GF2::_3d::InnerPrimitiveContainerT
                               , GF2::_3d::PrimitiveContainerT
                               , PclCloudT
                               >( argc, argv);
    }
    else if ( GF2::console::find_switch(argc,argv,"--3D") )
    {
        std::cout << "running planes" << std::endl;
        return ransacCli< GF2::PointContainerT
                       , GF2::_3d::InnerPrimitiveContainerT
                       , GF2::_3d::PrimitiveContainerT
                       , PclCloudT
                       >( argc, argv );
    }
    else
    {
        return ransacCli< GF2::PointContainerT
                       , GF2::_2d::InnerPrimitiveContainerT
                       , GF2::_2d::PrimitiveContainerT
                       , PclCloudT
                       >( argc, argv );
    } // if 3D

    return EXIT_FAILURE;
}

#undef TIC
#undef RETIC
#undef TOC
