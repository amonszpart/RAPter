#ifndef __GF2_POSTALIGNER_HPP__
#define __GF2_POSTALIGNER_HPP__

#include "post_opt/postAligner.h"
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include <limits>

namespace GF2
{
    template < class      PrimitivesT
               , class    PointsPtrT
               , typename Scalar
               , int      Dim >
    inline int
    PostAligner::run( PrimitivesT             & prims
                      , PointsPtrT              points
                      , std::vector<int> const& points_prims
                      , std::vector<int> const& prims_clusters_arg
                      , Eigen::Matrix<Scalar,3,1> (*getNormal)(typename PrimitivesT::value_type const& prim) )
    {
        using std::vector;
        using std::map;
        using std::set;

        map<int,int>  clusters_old_new;
        vector<int>   prims_clusters( prims_clusters_arg.size(), -1 ); // remapped

        for ( size_t prim_id = 0; prim_id != prims_clusters_arg.size(); ++prim_id )
        {
            const int cluster_id = prims_clusters_arg[prim_id];
            // insert new mapping
            if ( clusters_old_new.find(cluster_id) == clusters_old_new.end() )
                clusters_old_new[cluster_id] = clusters_old_new.size()-1;
            // apply mapping
            prims_clusters[ prim_id ] = clusters_old_new[cluster_id];

            std::cout << "cluster_id " << cluster_id << "->" << clusters_old_new[cluster_id] << std::endl;
        }
        const int C = clusters_old_new.size();

        // prims -> points
        vector< vector<int> > prims_points( prims.size() );
        for ( size_t pid = 0; pid != points->size(); ++pid )
        {
            prims_points[ points_prims[pid] ].push_back( pid );
        }

        // clusters -> primitives
        vector< vector<int> > clusters_primitives( C );
        for ( size_t prim_id = 0; prim_id != prims.size(); ++prim_id )
        {
            clusters_primitives[ prims_clusters[prim_id] ].push_back( prim_id );
        }

        //   p1.x p1.y p1.z 1
        // + p2.x
        Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> A( 1000, C * Dim + prims.size() );

        // prim_id -> x[j] (<nx,ny[,nz],d>)
        vector< vector<int> > xs;
        int curr_x = 0;

        // for each cluster
        for ( size_t cid = 0; cid != clusters_primitives.size(); ++cid )
        {
            Eigen::Matrix<Scalar,3,1> dir = getNormal( prims[clusters_primitives[cid][0]] );
            int curr_normal_start = curr_x;
            //  for each primitive
            for ( size_t prim_id_id = 0; prim_id_id != clusters_primitives[cid].size(); ++prim_id_id )
            {
                const int prim_id = clusters_primitives[cid][prim_id_id];

                // test
                if ( prims_clusters[prim_id] != cid ) std::cerr << "prims_clusters[prim_id] " << prims_clusters[prim_id] << " != " << cid << " cid" << std::endl;
                Eigen::Matrix<Scalar,3,1> dir2 = getNormal( prims[prim_id] );
                Scalar diff = (dir-dir2).norm();
                if  ( diff > 1e-3f )
                {
                    std::cerr << "diff[" << cid << "]: " << diff << " vs " << prim_id << "\t"
                              << dir.transpose() << " vs " << dir2.transpose() << std::endl;
                }

                // add mapping
                xs.push_back( vector<int>() );
                for ( int d = 0; d != Dim; ++d )
                {
                    xs.back().push_back( curr_normal_start + d );
                    if ( prim_id_id == 0 ) curr_x++;
                }
                xs.back().push_back( curr_x++ );

                //   add points
                for ( size_t pid_id = 0; pid_id != prims_points[prim_id].size(); ++pid_id )
                {
                    const int pid = prims_points[prim_id][pid_id];
                    // check
                    if ( points_prims[pid] != prim_id ) std::cerr << "points_prims[pid] " << points_prims[pid] << " != " << prim_id << " prim_id" << std::endl;

                    for ( int d = 0; d != Dim; ++d )
                    {
                        //A( xs[prim_id][d] )
                    }
                }
            }
        }

        for ( size_t i = 0; i != xs.size(); ++i )
        {
            for ( size_t d = 0; d != xs[i].size(); ++d )
            {
                std::cout << xs[i][d] << " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
        return EXIT_SUCCESS;
    } // ... run

    // get line inliers
    template< class PrimitivesT, class PointsPtrT
              , typename Scalar
              , int Dim        >
    int
    PostAligner::points2Prims( std::vector<int>                      & points_lines
                               , PointsPtrT                            cloud
                               , PrimitivesT                    const& prims
                               , float (*distance)(Eigen::Matrix<Scalar,3,1> const& pnt, Eigen::Matrix<Scalar,Dim,1> const& prim ) )

    {
        const int N = cloud->size();

        points_lines.resize( N, -1 );
        std::vector<double> min_distances( N, std::numeric_limits<float>::max() );

        // init sac model
        Scalar tmp_dist;

        // assignments - get closest line for each point
        for ( int prim_id = 0; prim_id != prims.size(); ++prim_id )
        {
            // calc distances
            for ( int pnt_id = 0; pnt_id != N; ++pnt_id )
            {
                if ( (tmp_dist = distance(cloud->at(pnt_id).getVector3fMap(), prims[prim_id]())) < min_distances[pnt_id] )
                {
                    min_distances[pnt_id] = tmp_dist;
                    points_lines [pnt_id] = prim_id;
                }
            }
        } // for line_id

        return EXIT_SUCCESS;
    } // get line inliers

} // ... ns GF2

#endif // __GF2_POSTALIGNER_HPP__
