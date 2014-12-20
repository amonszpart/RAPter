/**
 * FileName:    reassign.cpp
 * Author:      Aron Monszpart <a.monszpart@cs.ucl.ac.uk>
 * Created:     10/11/2014
 */

#define GCO_ENERGYTYPE float
#include "gco/GCoptimization.h"

namespace GF2
{
    template < typename _PointContainerT
             , class    _InnerPrimitiveContainerT
             , typename _PrimitiveContainerT
             , class    _PclCloudT >
    inline int reassignCli( int argc, char** argv )
    {
        typedef typename _PclCloudT::PointType                    PclPointT;
        typedef          std::map<int, _InnerPrimitiveContainerT> PrimitiveMapT;
        typedef typename _PointContainerT::value_type             PointPrimitiveT;
        typedef typename _InnerPrimitiveContainerT::value_type    PrimitiveT;
        typedef typename _PclCloudT::PointType                    PclPointT;
        typedef typename PrimitiveT::Scalar Scalar;

        std::cout << "hello reAssign\n";

        _PointContainerT         points;
        typename _PclCloudT::Ptr pcl_cloud;
        _PrimitiveContainerT     primitives;
        PrimitiveMapT            patches;
        RansacParams<Scalar>     params;

        // parse
        {
            bool valid_input = !parseInput<_InnerPrimitiveContainerT,_PclCloudT>( points, pcl_cloud, primitives, patches, params, argc, argv, false );

            if (     !valid_input
                  || (GF2::console::find_switch(argc,argv,"-h"    ))
                  || (GF2::console::find_switch(argc,argv,"--help")) )
            {
                std::cout << "[" << __func__ << "]: " << "Usage: " << argv[0] << "\n"
                          << "\t --cloud " << /*cloud_path <<*/ "\n"
                          << "\t -p,--prims " << /*input_prims_path <<*/ "\n"
                          << "\t -sc,--scale " << params.scale << "\n"
                          << "\t Example: ../ransac --assign --scale 0.03 --cloud cloud.ply -p patches.csv"
                          << "\n";

                return EXIT_FAILURE;
            }
        }

        int err = EXIT_SUCCESS;

        // debug visualize
        int plane_count = 0;
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
            vptr->addPointCloud<pcl::PointXYZ>( tmp.makeShared() );
            vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4.f );

            char name[255];
            for ( typename _PrimitiveContainerT::const_iterator it = primitives.begin(); it != primitives.end(); ++it )
            {
                for ( int lid = 0; lid != it->size(); ++lid )
                {
                    PrimitiveT const& prim = it->at( lid );
                    sprintf( name, "prim%06d_%06d", prim.getTag( PrimitiveT::GID ), lid );
                    vptr->addPlane( *(prim.modelCoefficients()), prim. pos()(0), prim. pos()(1), prim. pos()(2), name );

                    ++plane_count;
                }
            }

            vptr->spinOnce(150);
        }

        {
            const int num_pixels = pcl_cloud->size();
            const int num_labels = plane_count;

            Scalar *result = new Scalar[num_pixels];   // stores result of optimization

            // first set up the array for data costs
            std::map< int, std::pair<int,int> > labelMap;
            std::cout << "[" << __func__ << "]: " << "setting data labels..." << std::endl; fflush(stdout);
            Scalar *data = new Scalar[num_pixels*num_labels];

            #pragma omp for schedule(dynamic,8)
            for ( int pid = 0; pid < num_pixels; pid++ )
            {
                int label = 0;
                int lid0 = 0;
                for ( typename _PrimitiveContainerT::const_iterator it = primitives.begin(); it != primitives.end(); ++it, ++lid0 )
                    for ( int lid1 = 0; lid1 != it->size(); ++lid1, ++label )
                    {
                        Scalar dist = it->at(lid1).getDistance( points[pid] );
                        data[ pid*num_labels + label] = Scalar(100.) * dist * dist;
                        labelMap[ label ] = std::pair<int,int>( lid0, lid1 );
                    }
            }

            std::cout << "[" << __func__ << "]: " << "setting pairwise labels..." << std::endl; fflush(stdout);
            // next set up the array for smooth costs
            Scalar *smooth = new Scalar[ num_labels * num_labels ];
            #pragma omp for schedule(dynamic,8)
            for ( int l1 = 0; l1 < num_labels; l1++ )
                for (int l2 = 0; l2 < num_labels; l2++ )
                {
                    smooth[l1+l2*num_labels] = Scalar(1.);
                    //smooth[l1+l2*num_labels] = (l1-l2)*(l1-l2) <= 4  ? (l1-l2)*(l1-l2):4;
                }


            try
            {
                gco::GCoptimizationGeneralGraph *gc = new gco::GCoptimizationGeneralGraph(num_pixels,num_labels);
                gc->setDataCost(data);
                gc->setSmoothCost(smooth);

                // now set up a grid neighborhood system
                // first set up horizontal neighbors
                std::cout << "[" << __func__ << "]: " << "querying neighbourhood" << std::endl; fflush(stdout);
                std::vector< std::vector<int   > > neighs;
                std::vector< std::vector<Scalar> > sqr_dists;
                processing::getNeighbourhoodIndices( /*   [out] neighbours: */ neighs
                                                   , /* [in]  pointCloud: */ pcl_cloud
                                                   , /* [in]     indices: */ NULL
                                                   , /* [out]  sqr_dists: */ &sqr_dists
                                                   , /* [in]        nn_K: */ 15              // 15
                                                   , /* [in]      radius: */ params.scale    // 0.02f
                                                   , /* [in] soft_radius: */ false           // true
                                                   );

                std::cout << "[" << __func__ << "]: " << "setting neighbourhood" << std::endl; fflush(stdout);
                for ( int pid = 0; pid != neighs.size(); ++pid )
                    for ( int nid = 0; nid != neighs[pid].size(); ++nid )
                    {
                        Scalar d = Scalar(100.) * std::max( Scalar(0.), params.scale - sqr_dists[pid][nid]);
                        gc->setNeighbors( pid, neighs[pid][nid], d );
                    }


                printf("\nBefore optimization energy is %d",gc->compute_energy()); fflush(stdout);
                gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
                printf("\nAfter optimization energy is %d",gc->compute_energy());

                for ( int  i = 0; i < num_pixels; i++ )
                {
                    result[i] = gc->whatLabel(i);
                    std::pair<int,int> lidLid1 = labelMap[ result[i] ];
                    points[i].setTag( PointPrimitiveT::GID, primitives[lidLid1.first][lidLid1.second].getTag( PrimitiveT::GID) );
                }

                std::cout<<"result:";for(size_t vi=0;vi!=num_pixels;++vi)std::cout<<result[vi]<<" ";std::cout << "\n";
                std::cout << "[" << __func__ << "]: " << "writing to points_primitives.schnabel.csv" << std::endl;
                GF2::io::writeAssociations<PointPrimitiveT>( points, "./points_primitives.schnabel.csv" );

                delete gc;
            }
            catch (gco::GCException e)
            {
                e.Report();
            }

            delete [] result;
            delete [] smooth;
            delete [] data;
        }

    } //...reassignCli
} //...ns GF2
