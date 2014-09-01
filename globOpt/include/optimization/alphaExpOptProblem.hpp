#ifndef __GF2_ALPHAEXPOPTPROBLEM_HPP__
#define __GF2_ALPHAEXPOPTPROBLEM_HPP__

#include <vector>
#include <limits>

#include "localAnalysis.h"
#include "globFit2.h"


namespace am
{
    template <class PrimitivesT, class PointsPtrT>
    inline int AlphaExpOptProblem::solve( MaskType                           &  opt_mask
                                          , PrimitivesT                 const&  primitives
                                          , PointsPtrT                          cloud
                                          , Scalar                              working_scale
                                          , Eigen::Matrix<Scalar,-1,1>          lambdas
                                          , Scalar                      const   trunc_at
                                          , std::vector<Scalar>         const&  desired_angles
                                          , std::vector<int>                 *  labels_arg
                                          , Scalar                              angle_thresh
                                          )
    {
        typedef typename PointsPtrT::element_type::PointType PointT;
        //const Scalar gammasqr = 250.f;
        const int num_pixels = /*indices ? indices->size() :*/ cloud->size();
        const int num_labels = primitives.size();
        std::vector<int> labels;
        labels.resize( num_pixels );

        // init
        GCoptimizationGeneralGraph *gc = NULL;
        try
        {
            gc = new GCoptimizationGeneralGraph( num_pixels, num_labels );
            gc->setLabelCost ( lambdas(1) );
            std::cout << "set labelcost: " << lambdas(1) << std::endl;
        }
        catch ( GCException e )
        {
            e.Report();
        }

        /// DATA
        // first set up the array for data costs
        const float int_scale_sqr = 1e4f;
        int *data = new int[ num_pixels * num_labels ];
        {
            int capped_dists = 0, zero_dists = 0;
            float sum_zerod_dists = 0.f;
            for ( int pid = 0; pid != num_pixels; ++pid )
                for ( int line_id = 0; line_id != num_labels; ++line_id )
                {
                    float dist = primitives[line_id].point3Distance( cloud->at(pid).getVector3fMap() ) / working_scale;
                    float dist_original = dist;
                    dist *= lambdas(0) * dist * int_scale_sqr;

                    if ( (dist != dist) || (dist > 1e5) || (dist < 0) )
                    {
                        //std::cerr << "[" << __func__ << "]: " << "dist " << dist << std::endl;
                        ++capped_dists;
                        dist = 1e5;
                    }

                    //std::cout << "writing dist " << lambdas(0) * dist;
                    data[ pid * num_labels + line_id ] = dist;
                    if ( data[ pid * num_labels + line_id ] == 0 )
                    {
                        ++zero_dists;
                        sum_zerod_dists += dist_original;
                    }
                    //std::cout << " = " << data[ pid * num_labels + line_id ] << std::endl;
                }
            std::cout << std::endl;
            if ( capped_dists || zero_dists )
                std::cerr << "[" << __func__ << "]: "
                          << "capped dists: "     << capped_dists / static_cast<float>(num_pixels*num_labels) *100.f << "%, " << capped_dists << "/" << static_cast<float>(num_pixels*num_labels)
                          << ", zero dists:"      << zero_dists / static_cast<float>(num_pixels*num_labels)*100.f    << "%, " << zero_dists   << "/" << static_cast<float>(num_pixels*num_labels)
                          << ", lost distances: " << sum_zerod_dists<< ", avg: " << sum_zerod_dists / zero_dists
                          << std::endl;
            try
            {
                gc->setDataCost( data );
            }
            catch ( GCException e )
            {
                e.Report();
            }
        }

        /// SMOOTH
        // next set up the array for smooth costs
        const float cap = 1e5f;
        int *smooth = new int[ num_labels * num_labels ];
        {
            for ( int l1 = 0; l1 != num_labels; ++l1 )
                //for ( int l2 = 0; l2 != num_labels; ++l2 )
                for ( int l2 = l1; l2 != num_labels; ++l2 )
                {
                    Scalar angle = acos( primitives[l1].dir().dot(primitives[l2].dir()) );
                    // nan
                    if ( angle != angle )   angle = 0.f;

                    Scalar min_ang_diff = FLT_MAX, tmp;
                    for ( int ang_id = 0; ang_id != desired_angles.size(); ++ang_id )
                        if ( (tmp=fabs(angle - desired_angles[ang_id])) < min_ang_diff )
                            min_ang_diff = tmp;

                    Scalar ret = 0;
                    if ( min_ang_diff < trunc_at )
                    {
                        //ret = std::exp( 0.1f + 1.f * min_ang_diff ) - static_cast<Scalar>(1);
                        ret *= ret * ret;
                        ret *= ret;

                        //            ret = tmp * tmp * tmp;// * tmp * tmp * tmp
                        //            ret *= ret;
                    }

                    float val = lambdas(2) * ret * (l1 != l2) * int_scale_sqr * int_scale_sqr;
                    //std::cout << "setting smooth: " << val;
                    //if ( (val == 0.f) && (l1 != l2) ) val = 1;
                    if ( val > cap ) val = cap;
                    smooth[ l2 * num_labels + l1 ] = val;
                    smooth[ l1 * num_labels + l2 ] = val;

                    if ( smooth[ l2 * num_labels + l1 ] >= cap  )
                        std::cout << "smooth["<<l2<<","<<l1<<" ] = " << smooth[ l2 * num_labels + l1 ]
                                  << ", ang: " << min_ang_diff
                                  << ", ret: " << ret
                                  << ", val: " << val
                                  << std::endl;
                    //std::cout << " = " << smooth[ l2 * num_labels + l1 ] << std::endl;

                }

            try
            {
                gc->setSmoothCost( smooth );
            }
            catch ( GCException e )
            {
                e.Report();
            }

        }

#if 0
        // assign points to lines
        {
            using std::vector;
            vector< vector<int> > primitives_points;
            am::GlobFit2::assignPoints<PrimitivesT,PointsPtrT,Scalar>( /*        points_lines: */ static_cast< std::vector<int>* >(NULL)
                                                                       , /*              mask: */ static_cast< std::vector<int>* >(NULL)
                                                                       , /*        primitives: */ primitives
                                                                       , /*            points: */ cloud
                                                                       , /* primitives_points: */ &primitives_points
                                                                       , /*       dist_thresh: */ working_scale
                                                                       , /*         sqr_dists: */ static_cast<std::vector<Scalar>* >(NULL) );

            for ( size_t lid = 0; lid != primitives.size()-1; ++lid )
            {
                for ( size_t lid1 = lid+1; lid1 != primitives.size(); ++lid1 )
                {

                    std::pair<Scalar,Scalar> score__closest_angle = am::LocalAnalysis::matching( primitives[lid].dir(), primitives[lid1].dir(), desired_angles );
                    Scalar score = score__closest_angle.first;

                    if ( score < angle_thresh )
                    {
                        for ( size_t pid_id = 0; pid_id != primitives_points[lid].size(); ++pid_id )
                        {
                            const int pid = primitives_points[lid][pid_id];
                            for ( size_t pid_id2 = 0; pid_id2 != primitives_points[lid1].size(); ++pid_id2 )
                            {
                                const int pid2 = primitives_points[lid1][pid_id2];
                                try
                                {
                                    // (1-x^2)^2
                                    Scalar val = score * 1/angle_thresh;
                                    val *= val;
                                    val = (1-val);
                                    val *= val;
                                    //std::cout << "setting " << val << " for " << score << std::endl;
                                    gc->setNeighbors( std::min(pid,pid2), std::max(pid,pid2), val * int_scale_sqr );
                                }
                                catch ( GCException e )
                                {
                                    e.Report();
                                }
                            }
                        }
                    }
                }
                //std::cout << "lid" << lid << "/" << primitives.size() << "...finished" << std::endl;
            }
        }
#else
        std::vector<std::vector<int  > > neighs;
        std::vector<std::vector<float> > sqr_dists;
        smartgeometry::getNeighbourhoodIndices<PointT>( neighs
                                                        , cloud
                                                        , static_cast<std::vector<int> *>(NULL)
                                                        , &sqr_dists
                                                        , 10
                                                        , static_cast<float>( working_scale )
                                                        );

        // set neighbourhoods
        const Scalar working_scale_sqr = working_scale * working_scale;
        std::vector<int> neighvals; neighvals.reserve( neighs.size() * 15 );
        for ( size_t pid = 0; pid != cloud->size()-1; ++pid )
            for ( size_t pid2 = 1; pid2 != neighs[pid].size(); ++pid2 ) // don't count own
//             for ( size_t pid2 = pid+1; pid2 != cloud->size(); ++pid2 )
            {
                 float distsqr  = (cloud->at(pid).getVector3fMap() - cloud->at(pid2).getVector3fMap()).squaredNorm() * int_scale_sqr; //sqr_dists[pid][pid2] * int_scale_sqr;
//                    std::cout << "distsqr: " << distsqr
//                              << ", distsqr / working_scale_sqr : " << distsqr / working_scale_sqr
//                              << ", x = " << -1.f * distsqr / working_scale_sqr
//                              << ", exp( x ): " << exp( -1.f * distsqr / working_scale_sqr )
//                              << std::endl;
                 float neighval = distsqr / working_scale_sqr; //exp( -1.f * distsqr / working_scale_sqr ); // TODO: flip this to -1.f...
                 //std::cout << "setting neighval: " << neighval << " = " << static_cast<int>(neighval) << std::endl;
                 //                    gc->setNeighbors( pid, neighs[pid][pid2], neighval );
                 if ( neighval > 1e4f ) neighval = 1e4f;
                 try
                 {
                     gc->setNeighbors( pid, pid2, neighval );
                 }
                 catch ( GCException e )
                 {
                     e.Report();
                 }

                 neighvals.push_back( neighval );
            }

        // debug
        {
            int sum = std::accumulate( neighvals.begin(), neighvals.end(), 0 );
            std::cout << "neighval avg: " << sum / (float)neighvals.size();
            std::sort( neighvals.begin(), neighvals.end() );
            std::cout << ", median: " << neighvals[ neighvals.size() / 2 ];
            int n_min = neighvals[0];
            std::cout << ", neighval_min: " << n_min << ", neighval_max: " << neighvals.back() << std::endl;
        }

#endif
        try
        {
            //gc = new GCoptimizationGeneralGraph( num_pixels, num_labels );
            //gc->setDataCost  ( data       );
            //gc->setSmoothCost( smooth     );
            std::cout << "computing first energy...";
            printf("\tBefore optimization energy is %lld", gc->compute_energy() );
            gc->expansion();// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
            //gc->swap(5);
            printf("\tAfter optimization energy is %lld\n",gc->compute_energy());

            // copy output
            for ( int pid = 0; pid != num_pixels; ++pid )
            {
                labels[ pid ] = gc->whatLabel( pid );
            }

            // cleanup
            if ( gc     ) { delete gc; gc = NULL; }
        }
        catch ( GCException e )
        {
            e.Report();
        }
        catch ( std::exception &e )
        {
            std::cerr << "[" << __func__ << "]: " << e.what() << std::endl;
        }

        if ( labels_arg )
            *labels_arg = labels;

        if ( data   ) { delete [] data  ; data   = NULL; }
        if ( smooth ) { delete [] smooth; smooth = NULL; }

        return EXIT_SUCCESS;
    }

} // ns am

#endif // __GF2_ALPHAEXPOPTPROBLEM_HPP__
