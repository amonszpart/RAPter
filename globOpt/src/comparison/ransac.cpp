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
#endif // GF2_USE_PCL

template <typename _Scalar>
struct RansacParams
{
        enum ALGO { RANSAC };

    _Scalar scale;
    pcl::SacModel modelType;
    ALGO algo = RANSAC;
};

template < typename _PointContainerT
         , class    _InnerPrimitiveContainerT
         , typename _PrimitiveContainerT>
int ransacCli( int argc, char **argv )
{
    typedef typename _InnerPrimitiveContainerT::value_type    PrimitiveT;
    typedef typename _PointContainerT::value_type             PointPrimitiveT;
    typedef          std::map<int, _InnerPrimitiveContainerT> PrimitiveMapT;
    typedef          pcl::PointNormal                         PclPointT;
    typedef          pcl::PointCloud<PclPointT>               PclCloudT;
    typedef typename PointPrimitiveT::Scalar                  Scalar;

    bool valid_input = true;

    std::string cloud_path       = "./cloud.ply",
                input_prims_path = "./patches.csv",
                associations_path = "./points_primitives.csv";

    //am::Pearl::Params params;
    RansacParams<Scalar> params;
    std::string algo_str = "ransac";

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

        if (    (pcl::console::parse_argument( argc, argv, "-a", associations_path) < 0)
             && (pcl::console::parse_argument( argc, argv, "--assoc", associations_path) < 0)
             && (!boost::filesystem::exists(associations_path)) )
        {
            std::cerr << "[" << __func__ << "]: " << "-a or --assoc is compulsory" << std::endl;
            valid_input = false;
        }

        // scale
        if (    (GF2::console::parse_argument( argc, argv, "--scale", params.scale) < 0)
             && (GF2::console::parse_argument( argc, argv, "-sc"    , params.scale) < 0) )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale is compulsory" << std::endl;
            valid_input = false;
        }

        params.modelType = pcl::console::find_switch( argc, argv, "--3D" ) ? pcl::SacModel::SACMODEL_PLANE
                                                                           : pcl::SacModel::SACMODEL_LINE;

        // weights
        //pcl::console::parse_argument( argc, argv, "--unary", params.lambdas(0) );
        //pcl::console::parse_argument( argc, argv, "--cmp"  , params.beta );
        //valid_input &= pcl::console::parse_argument( argc, argv, "--pw"  , params.lambdas(2) ) >= 0;

        pcl::console::parse_argument( argc, argv, "--algo", algo_str );
        if ( !algo_str.compare("ransac") )
            params.algo = RansacParams<Scalar>::RANSAC;
        else
            std::cerr << "[" << __func__ << "]: " << "unrecognized ransac algorithm" << std::endl;

        if (     !valid_input
              || (GF2::console::find_switch(argc,argv,"-h"    ))
              || (GF2::console::find_switch(argc,argv,"--help")) )
        {
            std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << "\n"
                      << "\t --cloud " << cloud_path << "\n"
                      << "\t -p,--prims " << input_prims_path << "\n"
                      << "\t -a,--assoc " << associations_path << "\n"
                      << "\t -sc,--scale " << params.scale << "\n"
                      << "\t --3D\n"
                      //<< "\t --pw " << params.lambdas(2) << "\n"
                      //<< "\t --cmp " << params.beta << "\n"
                      << "\t Example: ../ransac --scale 0.03 --cloud cloud.ply -p patches.csv -a points_primitives.csv"
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
    PrimitiveMapT        patches, out_prims;
    if ( EXIT_SUCCESS == err )
    {
        std::cout << "[" << __func__ << "]: " << "reading primitives from " << input_prims_path << "...";
        GF2::io::readPrimitives<PrimitiveT, _InnerPrimitiveContainerT>( initial_primitives, input_prims_path, &patches );
        std::cout << "reading primitives ok (#: " << initial_primitives.size() << ")\n";
    } //...read primitives

    // read assoc
    std::vector<std::pair<int,int> > points_primitives;
    GF2::io::readAssociations( points_primitives, associations_path, NULL );
    for ( size_t i = 0; i != points.size(); ++i )
    {
        // store association in point
        points[i].setTag( PointPrimitiveT::GID, points_primitives[i].first );
    }


    std::vector<int> inliers;

    // created RandomSampleConsensus object and compute the appropriated model
    GF2::GidPidVectorMap populations;
    GF2::processing::getPopulations( populations, points );

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
                    pcl::SampleConsensusModelPlane<PclPointT>::Ptr model( new pcl::SampleConsensusModelPlane<PclPointT> (pcl_cloud, populations[gid]) );
                    pcl::RandomSampleConsensus<PclPointT> ransac (model);
                    ransac.setDistanceThreshold (params.scale);
                    if ( ransac.computeModel() )
                    {
                        ransac.getInliers(inliers);
                        // save line
                        if ( inliers.size() )
                        {
                            std::cout << "ret: " << ransac.model_coefficients_.transpose() << std::endl;

                            Eigen::Matrix<Scalar,3,1> nrm = ransac.model_coefficients_.head<3>();
                            nrm.normalize();
                            GF2::containers::add( out_prims, gid, PrimitiveT( Eigen::Matrix<Scalar,3,1>::Zero() - nrm * ransac.model_coefficients_(3), nrm) )
                                    .setTag( PrimitiveT::GID    , gid )
                                    .setTag( PrimitiveT::DIR_GID, gid );
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
                    std::cout << "[" << __func__ << "]: " << "doing ransac lines" << std::endl;
                    std::cout<<"populations["<<gid<<"]:";for(size_t vi=0;vi!=populations[gid].size();++vi)std::cout<<populations[gid][vi]<<" ";std::cout << "\n";
                    pcl::SampleConsensusModelLine<PclPointT>::Ptr model( new pcl::SampleConsensusModelLine<PclPointT> (pcl_cloud, populations[gid]) );

                    pcl::RandomSampleConsensus<PclPointT> ransac( model );
                    ransac.setDistanceThreshold( params.scale );
                    if ( ransac.computeModel() )
                    {
                        ransac.getInliers( inliers );

                        // save line
                        if ( inliers.size() )
                        {
                            GF2::containers::add( out_prims, gid, PrimitiveT( ransac.model_coefficients_.head<3>(), ransac.model_coefficients_.segment<3>(3)) )
                                    .setTag( PrimitiveT::GID    , gid )
                                    .setTag( PrimitiveT::DIR_GID, gid );
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
        }
    }


    GF2::io::savePrimitives<PrimitiveT,typename _InnerPrimitiveContainerT::const_iterator>( out_prims, "./primitives.ransac.csv" );
    GF2::io::writeAssociations<PointPrimitiveT>( points, "./points_primitives.ransac.csv" );


#if 0
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
                            , /* [in]        patches: */ &initial_primitives
                            );

        PrimitiveMapT out_prims;
        std::cout<<"labels:";for(size_t vi=0;vi!=labels.size();++vi)std::cout<<labels[vi]<<" ";std::cout << "\n";
        for ( int pid = 0; pid != primitives.size(); ++pid )
        {
            const int gid = labels[pid];
            // assign point
            points[pid].setTag( PointPrimitiveT::GID, gid );

            // add primitive
            if ( out_prims.find(gid) == out_prims.end() )
            {
                std::cout << "[" << __func__ << "]: " << "adding primitive[" << gid << "]: " << primitives[gid].toString() << std::endl;
                GF2::containers::add( out_prims, gid, primitives[gid] )
                        .setTag( PrimitiveT::GID    , gid )
                        .setTag( PrimitiveT::DIR_GID, gid );

            }
        }

        GF2::io::writeAssociations<PointPrimitiveT>( points, "./points_primitives.pearl.csv" );
        GF2::io::savePrimitives<PrimitiveT,typename _InnerPrimitiveContainerT::const_iterator>( out_prims, "./primitives.pearl.csv" );
    } //...work
#endif

    return err;
}


int main(int argc, char *argv[])
{
    std::cout << "hello ransac\n";
    if ( GF2::console::find_switch(argc,argv,"--schnabel3D") )
    {
        typedef GF2::PointPrimitive PointPrimitiveT;
        typedef GF2::PlanePrimitive PrimitiveT;
        typedef pcl::PointNormal    PT;

        std::string cloud_path = "./cloud.ply";
        std::cout << "hello Schnabel\n";
        if ( (  GF2::console::parse_argument( argc, argv, "--cloud", cloud_path) < 0)
             && (!boost::filesystem::exists(cloud_path)) )
        {
            std::cerr << "--cloud" << std::endl;
            return 1;
        }

        // scale
        float scale = 0.1f;
        if (    (GF2::console::parse_argument( argc, argv, "--scale", scale) < 0)
             && (GF2::console::parse_argument( argc, argv, "-sc"    , scale) < 0) )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale is compulsory" << std::endl;
            return 1;
        }

        int min_support_arg = 300;
        if (    (GF2::console::parse_argument( argc, argv, "--minsup", min_support_arg) < 0) )
        {
            std::cerr << "[" << __func__ << "]: " << "--minsup required (300?)" << std::endl;
            return EXIT_FAILURE;
        }

        int err = EXIT_SUCCESS;

        std::vector<PointPrimitiveT> points;
        pcl::PointCloud<PT>::Ptr   pcl_cloud( new pcl::PointCloud<PT> );
        if ( EXIT_SUCCESS == err )
        {
            err = GF2::io::readPoints<PointPrimitiveT>( points, cloud_path, &pcl_cloud );
            if ( err != EXIT_SUCCESS )  std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
        } //...read points
        pcl::PointCloud<pcl::PointXYZ> tmp;
        for ( size_t pid = 0; pid != pcl_cloud->size(); ++pid )
        {
            tmp.push_back( pcl::PointXYZ( pcl_cloud->at(pid).x
                                          , pcl_cloud->at(pid).y
                                          , pcl_cloud->at(pid).z) );
        }
        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
        vptr->setBackgroundColor( .5, .6, .6 );
        vptr->addPointCloud<pcl::PointXYZ>( tmp.makeShared(), "asdf");
        vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4.f );

        vptr->spin();

        for ( size_t i = 0; i != std::min(10UL,points.size()); ++i )
        {
            std::cout << "\tpoints[" << i << "]: " << points[i].toString() << "\n\tcloud[" << i << "]: " << pcl_cloud->at(i).getVector3fMap().transpose() << std::endl;
        }

         std::vector<GF2::PlanePrimitive> planes;
         typedef std::map<int,int> PidGidT;
         PidGidT pidGid;

        err = GF2::SchnabelEnv::run<pcl::PointCloud<PT> >( planes
                                   , pidGid
                                   , points
                                   , pcl_cloud
                                   , scale
                                   , min_support_arg  );
        std::map< int, std::vector<GF2::PlanePrimitive> > out_prims;
        for ( int gid = 0; gid != planes.size(); ++gid )
        {
            GF2::containers::add( out_prims, gid, planes[gid] )
                    .setTag( PrimitiveT::GID    , gid )
                    .setTag( PrimitiveT::DIR_GID, gid );
            std::cout << "added " << out_prims[gid].back().toString() << std::endl;
        }
#if 0
        // assign points
        std::cout << "starting assignment" << std::endl; fflush( stdout );
        for ( int pid = 0; pid != points.size(); ++pid )
        {
            float min_dist = FLT_MAX, tmp; int min_gid = 0;
            for ( int gid = 0; gid != planes.size(); ++gid )
            {
                tmp = planes[gid].getDistance( points[pid].pos() );
                if ( (tmp < min_dist) && (tmp < scale) )
                {
                    min_dist = tmp;
                    min_gid = gid;
                }
            }
            pidGid[ pid ] = min_gid;
            points[pid].setTag( PointPrimitiveT::GID, min_gid );
        }
        std::cout << "finishing assignment" << std::endl;
#endif
        // debug assignments
        pcl::PointXYZ pnt, plane_pnt;
        for ( int pid = 0; pid != points.size(); ++pid )
        {
            if ( !(pid % 1000) )
            {
                char name[255];
                sprintf( name, "pointarrow%06d", pid );
                pnt.getVector3fMap() = points[pid].pos();
                const int gid = points[pid].getTag( PointPrimitiveT::GID );
                plane_pnt.getVector3fMap() = out_prims[gid].at(0).pos();
                vptr->addArrow( plane_pnt, pnt,  .1, .9, .5, false, name );
            }
        }

//        for ( PidGidT::const_iterator it = pidGid.begin(); it != pidGid.end(); ++it )
//        {
//            points[it->first].setTag( PointPrimitiveT::GID, it->second );
//        }
        for ( size_t i = 0; i != std::min(6UL,planes.size()); ++i )
        {
            char name[255];
            sprintf( name, "plane%d", i );
            vptr->addPlane( *(planes[i].modelCoefficients()), planes[i]. pos()(0), planes[i]. pos()(1), planes[i]. pos()(2), name );
        }
        vptr->addCoordinateSystem(0.5);
        vptr->spin();

        std::cout << " writing " << out_prims.size() << " primitives to ./primitives.schnabel.csv" << std::endl;
        GF2::io::savePrimitives<PrimitiveT,std::vector<GF2::PlanePrimitive>::const_iterator >( out_prims, "./primitives.schnabel.csv" );
        GF2::io::writeAssociations<PointPrimitiveT>( points, "./points_primitives.schnabel.csv" );

        return EXIT_SUCCESS;


    }
    else if ( GF2::console::find_switch(argc,argv,"--3D") )
    {
        std::cout << "running planes" << std::endl;
        return ransacCli< GF2::PointContainerT
                       , GF2::_3d::InnerPrimitiveContainerT
                       , GF2::_3d::PrimitiveContainerT
                       >( argc, argv );
        return EXIT_FAILURE;
    }
    else
    {
        return ransacCli< GF2::PointContainerT
                       , GF2::_2d::InnerPrimitiveContainerT
                       , GF2::_2d::PrimitiveContainerT
                       >( argc, argv );
    } // if 3D

    return EXIT_FAILURE;
}

