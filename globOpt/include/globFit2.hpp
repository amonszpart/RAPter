#ifndef __GF2_GLOBFIT2_HPP__
#define __GF2_GLOBFIT2_HPP__

#include <chrono>
#include <limits>
#include "optimization/simAnnOptProblem.h"
#include "optimization/alphaExpOptProblem.h"
#include "AMUtil2.h"
#include "pcltools/util.hpp"
#include "localAnalysis.h"
#include "pcl/search/kdtree.h"

#define TIC auto start_timing = std::chrono::system_clock::now();
#define RETIC start_timing = std::chrono::system_clock::now();
#define TOC(title) { std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_timing; \
                     std::cout << title << ": " << elapsed_seconds.count() << " s" << std::endl; }
#define CLUSTERED_INIT 1

template <typename PointT, typename Scalar> inline PointT
pntToPCL( Eigen::Matrix<Scalar,3,1> const& pnt )
{
    PointT ret;
    ret.x = pnt(0); ret.y = pnt(1); ret.z = pnt(2);

    return ret;
}

template <class PrimitivesT> int
am::GlobFit2::optimize( MaskType                                  & min_config
                        , PrimitivesT                        const& lines
                        , PrimitiveClustering<PrimitivesT>   const* clustering
                        , pcl::PointCloud<MyPoint>::Ptr             cloud
                        , Eigen::Matrix<float,-1,1>          const& lambdas
                        , float                              const  scale
                        , std::vector<float>                 const& desired_angles
                        , int                                const  max_iterations
                        , int                                const  max_step_count
                        , double                             const  threshold       // line inlier threshold
                        , float                              const  trunc_pw_at_angle // in radians
                        , std::vector<int>                   const* indices
                        , int                                const  fixedK
                        , std::vector<int>                        * labels
                        , float                                   * min_e_arg
                        , std::string                        const* p_out_dir )
{
    typedef typename PrimitivesT::value_type::Scalar Scalar;

    using std::vector;

    if ( !clustering ) { std::cerr << "[" << __func__ << "]: " << " no clustering...can't work\n" << std::endl; return EXIT_FAILURE; }

    std::cout << "[" << __func__ << "]: " << "start"
              << ", cloud->size(): " << cloud->size()
              << ", lines.size(): " << lines.size()
              << ", clusters.size(): " << clustering->clusters_lines.size()
              << std::endl;

    if ( fixedK < 1 )
    {
        std::cerr << "[" << __func__ << "]: " << "fixedK < 1...exiting";
        return EXIT_FAILURE;
    }

    int                 min_step_id = 0; // when did it find the current solution?
    float               min_e       = FLT_MAX;
    std::vector<int>    min_step_ids( max_iterations, 0);
    std::vector<float>  min_es( max_iterations, FLT_MAX );
    std::vector<MaskType> min_configs( max_iterations, std::vector<int>(lines.size(),0) );
    std::vector<Eigen::Matrix<float,-1,1>> min_energies( max_iterations, Eigen::Vector4f::Zero() );

    std::vector<int> initial_config( clustering->clusters_lines.size(), 0 );
    int target_count = ceil(fixedK / static_cast<float>(clustering->directions_in_cluster) );
    int start        = std::max(0,(int)(initial_config.size()-target_count)/2);
    int stop         = std::min(initial_config.size()-1,(initial_config.size()+target_count)/2);
    if ( stop - start != fixedK ) { std::cerr << "correcting " << stop - start << " by " << fixedK - (stop-start) << "\n"; start -= fixedK - (stop-start); }
    start = std::max( start, 0 ); stop = std::min( stop, (int)initial_config.size()-1 );

    for ( int i = start; i != stop; ++i )   { if ( i >= initial_config.size() ) std::cerr << "werwer" << std::endl; initial_config[ i ] = 1; };

    if ( (std::accumulate( initial_config.begin(), initial_config.end(), 0) != fixedK) && (target_count < clustering->clusters_lines.size()) )
    {
        std::cerr << "[" << __func__ << "]: " << "initial config not correct number of clusters: " << fixedK
                  << " vs. " << std::accumulate( initial_config.begin(), initial_config.end(), 0 )
                  << std::endl;
    }

    // Energy stats
    vector<vector<float> > energies( max_iterations );
    for ( int i = 0; i != max_iterations; ++i ) energies[i].reserve( max_step_count );

    int finished_it_count = 0;

#if 0
    AlphaExpOptProblem::solve( min_config
                               , lines
                               , cloud
                               , scale
                               , lambdas
                               , trunc_pw_at_angle
                               , desired_angles
                               , labels );
    return EXIT_SUCCESS;

#endif

    TIC;
#   pragma omp parallel for num_threads(8) shared(min_e, min_config, min_step_id, energies)
    for ( int it = 0; it < max_iterations; ++it )
    {
        PointPrimitiveDistanceCache<Scalar> cache( lines.size(), 200 );
        SimAnnOptProblem<PrimitivesT> sa_problem( lines, cloud, threshold, indices, lambdas, scale, max_step_count, trunc_pw_at_angle, &cache );
        sa_problem._omp_thread_id = it;

        // INIT
        {
            int init_retries = 0;
            do
            {
                std::random_shuffle( initial_config.begin(), initial_config.end() );
            } while ( (EXIT_SUCCESS != sa_problem.initWithClusters(*clustering, initial_config, fixedK)) && (init_retries++ < 50) );

            if ( init_retries >= 50 )
            {
                std::cerr << "[" << __func__ << "][" << it << "]: "
                          << "sa_problem.initWithClusters failed 50 times\n";
                continue;
            }
        }

        // init
        sa_problem.getEnergy( &desired_angles );

        // run
        float e;
        while ( sa_problem.hasNext() )
        {
            // step
            sa_problem.next( &desired_angles );
            // read
            e = sa_problem.getEnergy( &desired_angles );
            // minimize
            if ( e < min_es[it] )
            {
                min_es[it]       = e;
                min_configs[it]  = sa_problem.getConfig();
                min_step_ids[it] = sa_problem.getStepCount();
            }

            // stats
            energies[it].push_back( e );
        }

        // print minimum THREAD energy
        {
            sa_problem.calcEnergy( min_configs[it], desired_angles, &(min_energies[it]) );
        } // ... print

        // statistics
        {
            sa_problem.report();
        }

#       pragma omp critical
        {
            std::cout << "[" << __func__ << "]: " << "progress: " << ++finished_it_count << "/" << max_iterations << std::endl;
        }
    } // ... for iterations

    // gather over threads
    {
        int min_pos = std::distance( min_es.begin(),std::min_element(min_es.begin(), min_es.end()) );
        min_e       = min_es      [ min_pos ];
        min_config  = min_configs [ min_pos ];
        min_step_id = min_step_ids[ min_pos ];

        // print minimum OVERALL energy
        {
            std::cout << "MIN_CONFIG: " << min_energies[min_pos].transpose()
                      << ", found it at iteration: " << min_step_id
                      << ", energy: " << min_e << " = ";
            for ( int lambda_id = 0; lambda_id != lambdas.rows(); ++lambda_id )
                std::cout << min_energies[min_pos](lambda_id) * lambdas(lambda_id) << ((lambda_id==lambdas.rows()-1)?"":" + ");
            std::cout << std::endl;
        } // ... print energy

        if ( min_e_arg ) *min_e_arg = min_e;
    } // ... gather


    TOC(__func__);

    // gather stats
    {
        std::ofstream f( *p_out_dir + "energies.txt" );
        for ( size_t step = 0; step != max_step_count; ++step )
        {
            float min_energy = FLT_MAX;
            for ( size_t it = 0; it != energies.size(); ++it )
            {
                f << it << " " << step << " " << energies[it][step] << "\n";
                if ( energies[it][step] < min_energy ) min_energy = energies[it][step];
            }
            f << energies.size() /* == last_it+1 */ << " " << step << " " << min_energy << std::endl;
        }

        f.close();

        // histogram
        {
            std::vector<unsigned> hist(1000,0);
            for ( size_t step = 0; step != max_step_count; ++step )
            {
                for ( size_t it = 0; it != energies.size(); ++it )
                {
                    unsigned h = (energies[it][step] - fixedK) * 10000.f;
                    if ( h >= hist.size() ) ++hist.back();
                    else ++hist[h];
                }
            }
            std::ofstream f_hist( *p_out_dir + "hist.txt" );
            for ( size_t hi = hist.size(); hi; --hi )
                f_hist << hist[hi-1] << std::endl;
            f_hist.close();
        }
    }

    // calculate labels as output
    if ( labels )
    {
        SimAnnOptProblem<PrimitivesT>::calcDataTerms( cloud, lines, min_config, indices, &scale, labels );
    }

    // print output
    {
        std::cout << "selection: ";
        for ( size_t line_id = 0; line_id != lines.size(); ++line_id )
            if ( min_config[line_id] )
                std::cout << line_id << ", ";
        std::cout << std::endl;
    }

    std::cout << "[" << __func__ << "]: " << "finished " << fixedK << std::endl; fflush(stdout);
    return EXIT_SUCCESS;
} // ... optimize

// assigns each point to a line
// also each line points inside scale
template <class PrimitivesT, class PointsPtrT, typename Scalar> inline int
am::GlobFit2::assignPoints( std::vector<int>          * points_lines
                            , MaskType           const* mask
                            //, GF2::SelectionType const& selection
                            , PrimitivesT        const& primitives
                            , PointsPtrT                cloud
                            , std::vector< std::vector<int> > *primitives_points
                            , Scalar             const  dist_threshold // to assign points to primitives
                            , std::vector<Scalar>     * p_distances_sqr )
{
    const int N = cloud->size();
    std::vector<Scalar> min_distances_sqr( N, std::numeric_limits<Scalar>::max() );
    if ( points_lines ) points_lines->resize( N, 0 );

    if ( primitives_points && (primitives.size() > primitives_points->size()) )
        primitives_points->resize( primitives.size() );

    Scalar tmp_dist_sqr, dist_threshold_sqr = dist_threshold * dist_threshold;
    for ( size_t pnt_id = 0; pnt_id != N; ++pnt_id )
    {
        //GF2::SelectionType::const_iterator end_it = selection.end();
        //for ( GF2::SelectionType::const_iterator it = selection.begin(); it != end_it; ++it )
        typename PrimitivesT::const_iterator end_it = primitives.end();
        int prim_id = 0;
        for ( typename PrimitivesT::const_iterator it = primitives.begin(); it != end_it; ++it, ++prim_id )
        {
            if ( mask && !((*mask)[prim_id]) ) continue;

            tmp_dist_sqr  = it->point3Distance( cloud->at(pnt_id).getVector3fMap() );
            tmp_dist_sqr *= tmp_dist_sqr;

            if ( points_lines && tmp_dist_sqr < min_distances_sqr[pnt_id] )
            {
                min_distances_sqr[pnt_id] = tmp_dist_sqr;
                (*points_lines)  [pnt_id] = prim_id;
            }

            // assign point to lines, if inside scale
            if ( primitives_points && (tmp_dist_sqr < dist_threshold_sqr) )
            {
                (*primitives_points)[ prim_id ].push_back( pnt_id );
            }
        }
    }

    if ( p_distances_sqr )
        *p_distances_sqr = min_distances_sqr;

    return EXIT_SUCCESS;
}

template <class PointsPtrT, typename Scalar = float, int Dim> inline int
getDistancesOnLine( std::vector<Scalar>                    & distances
                    , Eigen::Matrix<Scalar,Dim,1>     const& line
                    , PointsPtrT                             cloud
                    , Eigen::Matrix<Scalar,2,1>            * minMax = NULL )
{
    // calculate projected points' distances from P0 along the normal
    Scalar min_val = DBL_MAX, max_val = -DBL_MAX;  // farthest point distance along Normal "behind" P0, farthest point distance along Normal from P0

    distances.resize( cloud->size() );
    // all projected points
    for ( int point_id = 0; point_id < cloud->size(); ++point_id )
    {
        // get vector P1-P0
        Eigen::Matrix<Scalar,3,1> diff = (cloud->at(point_id).getVector3fMap() - line.template head<3>() );
        // length
        distances[point_id] = diff.norm();
        // direction
        if ( diff.dot( line.template head<3>() + line.template segment<3>(3)) < 0.f ) // if dot product doesn't align with local normal
            distances[point_id] *= -1.f;                            // distance is behind (negative)
        // min, max
        if ( distances[point_id] > max_val )    max_val = distances[point_id];
        if ( distances[point_id] < min_val )    min_val = distances[point_id];
    }

    if ( minMax )
    {
        (*minMax)(0) = min_val;
        (*minMax)(1) = max_val;
    }

    return EXIT_SUCCESS;
}

template <class PointsPtrT, typename Scalar = float, int Dim> inline int
lineHistogram( std::vector<Eigen::Matrix<Scalar,3,1> > & ends
               , Eigen::Matrix<Scalar,Dim,1>             line
               , PointsPtrT                              cloud
               , int                                     N_bins    = 10
               , int                                     min_limit = 5
             )
{
    std::vector<Scalar> distances;
    Eigen::Matrix<Scalar,2,1> minMax;
    getDistancesOnLine( distances, line, cloud, &minMax );
    Scalar min_val = minMax(0), max_val = minMax(1);

    // histogram bin width

    int beg_pos = 0;
    int end_pos = N_bins-1;
    double step;
    std::vector<int> histogram;

    bool needRestart     = true;
    int  restarted_count = -1;
    while ( needRestart && (restarted_count < 1) )
    {
        needRestart = false;
        ++restarted_count;

        step = (max_val - min_val) / N_bins;
        // prepare histogram
        histogram.resize( N_bins, 0 );

        // bin points
        for ( int point_id = 0; point_id < cloud->size(); ++point_id )
        {
            // calculate bin index
            int index = (int) floor((distances[point_id] - min_val) / step);
            // check, if histogram has been clamped, and this doesn't fit
            if ( (index < 0) || (index >= histogram.size()) )
                continue;

            // increment bin
            ++histogram[ index ];
        }

        int median = am::util::algo::median( histogram );

        // check histogram, recurr if necessary
        if ( (beg_pos == 0) && (end_pos == histogram.size()-1) )
        {
            // calculate minimum bin size limit
            int limit;
            {
                limit = std::max( min_limit, median/2 );
                std::cout << "limit is: " << limit << std::endl;
            }

            // check from beginning for first >limit
            while ( ( beg_pos < histogram.size()          ) && (histogram[beg_pos] < limit) )      ++beg_pos;
            // check from end for first <limit
            while ( ((end_pos > 0) && (beg_pos < end_pos) ) && (histogram[end_pos] < limit) )      --end_pos;

            // if there is something to clamp, restart
            if ( (beg_pos > 0) || (end_pos < histogram.size() - 1) )
            {
                if ( beg_pos >= end_pos )
                {
                    std::cerr << "not rescaling, since beg>= end" << std::endl;
                    break;
                }

                max_val = min_val + step * (end_pos + 1);
                min_val += step * beg_pos; // histogram.mins[beg_pos];
                std::cerr << "new min_val: " << min_val << ", max_val: " << max_val << std::endl;

                needRestart = true;
            }
        }
    }
    ends.push_back( (Eigen::Matrix<Scalar,3,1>() << min_val + beg_pos * step, 0, 0).finished() );
    ends.push_back( (Eigen::Matrix<Scalar,3,1>() << min_val + (end_pos+1) * step, 0, 0).finished() );

    return EXIT_SUCCESS;
}

/**
 * @brief getDistanceIntervalsSorted pair<Scalar,Scalar>: pair.first is the interval length between pair.second and the one that comes after pair.second in sorted_distances
 * @param distanceIntervals
 * @param distances
 */
template <class PointsPtrT, typename Scalar> inline int
getDistanceIntervalsSorted( std::vector<std::pair<Scalar,Scalar>>  & distanceIntervals
                            , std::vector<Scalar>                  & distances // need the copy to sort
                            , PointsPtrT                 /*cloud*/              )
{
    using std::vector;
    using std::pair;
    std::sort( distances.begin(), distances.end() );
    distanceIntervals.resize( distances.size() - 1 );
    for ( size_t did = 0; did != distances.size()-1; ++did )
    {
        distanceIntervals[did] = std::pair<Scalar,Scalar>( distances[did+1] - distances[did], distances[did] );
    }

    return EXIT_SUCCESS;
}

struct Gap
{
        Gap() : start(0), end(0), length(0) {}
        Gap( int start, int end, int length ) : start(start), end(end), length(length) {}
        int start = 0, end = 0, length = 0;
};

template <typename Scalar> inline int
getLongestGap( std::vector<Scalar>                           & gapBegEnd
               , std::vector<std::pair<Scalar,Scalar> > const& distanceIntervals // ! distanceIntervals has one less members than cloud
               , Scalar                                 const  gap_threshold )
{
    Gap best;
    std::vector<Gap> gaps;

    int gap_start = 0, gap_end = 0;
    for ( size_t did = 0; did != distanceIntervals.size(); ++did )
    {
        if ( distanceIntervals[did].first < gap_threshold )
            gap_end = did; // extend current gap
        else
        {
            // check if long enough
            if ( gap_end - gap_start > best.length )
            {
                best.length = gap_end - gap_start;
                best.start  = gap_start;
                best.end    = gap_end;
            }

            gaps.push_back( (Gap){gap_start, gap_end, gap_end - gap_start} );

            // reset
            gap_start = gap_end = did + 1;
        }
    }

    // check, if best was written
    if ( gap_end - gap_start > best.length )
    {
        best.length = gap_end - gap_start;
        best.start  = gap_start;
        best.end    = gap_end;
    }

    // output
    if ( best.start >= best.end )
        std::cerr << "[" << __func__ << "]: " << "wtf...start >= end: " << best.start << ", " << best.end << std::endl;
    else
    {
        //gapBegEnd = { distanceIntervals[best.start].second, distanceIntervals[best.end].second + distanceIntervals[best.end].first  };
        for ( int i = 0; i != gaps.size(); ++i )
        {
            if ( gaps[i].length <= 1 ) continue;

            gapBegEnd.push_back( distanceIntervals[gaps[i].start].second );
            gapBegEnd.push_back( distanceIntervals[gaps[i].end  ].second + distanceIntervals[gaps[i].end].first );
        }
    }
    return EXIT_SUCCESS;
}

template <class PrimitivesT, class PointsPtrT, typename Scalar> inline int
am::GlobFit2::split( std::vector<std::vector<Scalar > > &lines_ends
                     , am::MaskType const& mask
                     , PrimitivesT  const& prims
                     , PointsPtrT          cloud
                     , Scalar       const  scale
                     , std::vector<int>  * points_lines_arg )
{
    typedef typename PointsPtrT::element_type::PointType PointT;
    typedef typename PrimitivesT::value_type PrimitiveT;

    // convert to sparse set
    GF2::SelectionType selection;
    {
        for ( int prim_id = 0; prim_id != prims.size(); ++prim_id )
        {
            if ( !mask[prim_id] ) continue;
            selection.insert( prim_id );
        }
    }
    if ( !selection.size() )
    {
        std::cerr << "[" << __func__ << "]: " << "no candidates choosen..." << std::endl;
        return EXIT_FAILURE;
    }

    // assign points
    std::vector<int> points_lines;
    {
        assignPoints( &points_lines, &mask, prims, cloud );
    }

    // project
    pcl::SampleConsensusModel<PointT> *sac_model = NULL;
    if ( PrimitiveT::Dim == 6 )     sac_model = new pcl::SampleConsensusModelLine <PointT>( cloud );
    else                            sac_model = new pcl::SampleConsensusModelPlane<PointT>( cloud );

    //        std::vector<int> neighbour_indices;
    //        std::vector<float> sqr_dists;

    // prepare output
    lines_ends.resize( selection.size() );

    // for every primitive
    int prim_id = 0;
    GF2::SelectionType::const_iterator end_it = selection.end();
    for ( GF2::SelectionType::const_iterator it = selection.begin(); it != end_it; ++it, ++prim_id )
    {
        // collect inlier points
        std::vector<int> indices; indices.reserve( cloud->size() / selection.size() );
        for ( size_t point_id = 0; point_id != points_lines.size(); ++point_id )
        {
            if ( *it == points_lines[point_id] )
                indices.push_back( point_id );
        }
        if ( indices.size() < 2 ) continue;

        // project points
        PointsPtrT projected_cloud( new typename PointsPtrT::element_type() );
        std::cout << "prims.size(): "  << prims.size() << ", *it: " << *it  << std::endl;
        std::cout << "prims[*it].coeffs():";for(size_t vi=0;vi!=prims[*it].coeffs().size();++vi)std::cout<<prims[*it].coeffs()[vi]<<" ";std::cout << "\n";
        std::cout << "projected_cloud->size: " << projected_cloud->size() << ", indices.size: " << indices.size() << ", indices[0]: " << indices[0] << std::endl;
        sac_model->projectPoints( indices, prims[*it].coeffs(), *projected_cloud, false );

        // histogram - deprecated
        //lineHistogram( lines_ends[prim_id], prims[*it](), projected_cloud, 10, 5 );

        // distances to P0
        std::vector<Scalar> distances;
        getDistancesOnLine( distances, prims[*it](), projected_cloud );

        // sort intervals
        std::vector<std::pair<Scalar,Scalar> > distanceIntervals; // contains an interval and a starting distance for it
        getDistanceIntervalsSorted( distanceIntervals, distances, projected_cloud );

        // find longest innterruptionless gap using threshold
        std::vector<Scalar> gapBegEnd;
        getLongestGap( gapBegEnd, distanceIntervals, scale );

        if ( gapBegEnd.size() )
            lines_ends[prim_id] = gapBegEnd;
        else
        {
            std::cout << "brace youreself";
            lines_ends[prim_id] = { distanceIntervals.front().second, distanceIntervals.back().second + distanceIntervals.back().first };
            std::cout << " huh" << std::endl;
        }
    }

    if ( sac_model ) { delete sac_model; sac_model = NULL; }

    if ( points_lines_arg )
        *points_lines_arg = points_lines;

    return EXIT_SUCCESS;
}

//template <class PrimitivesT, class PointsPtrT> inline int
//am::GlobFit2::second( PrimitivesT         & out
//                     , am::MaskType const& mask
//                     , PrimitivesT  const& prims
//                     , PointsPtrT          cloud )
//{

//}

#if 0
template <class PrimitivesT, typename Scalar> inline int
am::GlobFit2::labelPoints( MaskType                                  & opt_mask
                           , PrimitivesT                        const& primitives
                           , PrimitiveClustering<PrimitivesT>   const* clustering
                           , pcl::PointCloud<MyPoint>::Ptr             cloud
                           , float                              const  scale      )
{

}
#endif

template <class PrimitivesT> int
am::GlobFit2::visualize( PrimitivesT                          const& lines
                         , MaskType                           const& optimized_mask
                         , pcl::PointCloud<MyPoint>::Ptr             cloud
                         , float                              const  scale
                         , std::vector<float>                 const& desired_angles
                         , MaskType                           const* real_mask
                         , pcl::PointIndices::Ptr                    indices
                         , std::vector<int>                   const* labels
                         , bool                               const  spin
                         , bool                               const  show_angles
                         , float                              const  gap_threshold )
{
    typedef typename PrimitivesT::value_type PrimitiveT;
    typedef typename PrimitiveT::Scalar Scalar;
    const bool show_line_id = false;

    int dim = PrimitiveT::Dim;
    const int NCorners = ( (dim==6) ? 2 : 4 );

    using std::vector;
    std::cout << "[" << __func__ << "]: " << "lines.size(): " << lines.size() << std::endl;
    std::cout << "[" << __func__ << "]: " << "optimized_mask:" << std::accumulate( optimized_mask.begin(), optimized_mask.end(), 0) << std::endl;

    auto getCloudName    = [](             int vpid) { return am::util::sprintf("cloud%d",vpid); };
    auto getLineName     = [](int line_id, int vpid) { return am::util::sprintf<int>("line%d_vp%d",line_id,vpid); };
    auto getLineTextName = [](int line_id, int vpid) { return am::util::sprintf<int>("line%d_tag_vp%d",line_id,vpid); };

    // init vis
    pcl::visualization::PCLVisualizer::Ptr vptrs[2] = { pcl::visualization::PCLVisualizer::Ptr(new pcl::visualization::PCLVisualizer()),
                                                        pcl::visualization::PCLVisualizer::Ptr(new pcl::visualization::PCLVisualizer()) };
    int vpids[2] = {0,0};
    // show cloud
    {
//            vptr->createViewPort( 0., 0., .5, 1., vpids[0] );
//            vptr->createViewPort( .5, 0., 1., 1., vpids[1] );
        for ( int vpid = 0; vpid != 2; ++vpid )
        {
            std::string cloud_name = getCloudName( vpid );
            vptrs[vpid]->setBackgroundColor( 1., 1., 1., vpids[vpid] );
            vptrs[vpid]->addPointCloud( cloud, cloud_name, vpids[vpid] );
            vptrs[vpid]->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3., cloud_name, vpids[vpid] );
            vptrs[vpid]->setPosition( vpid * 640, 30 );
            vptrs[vpid]->resetCamera();
        }
    }

    // show lines
    std::vector<pcl::PointXYZ> endpoints; endpoints.reserve( lines.size()*2 );

    // labels_points
    vector< vector<int> > labels_points( lines.size() );
    if ( labels )
        for ( size_t pid = 0; pid != (indices ? indices->indices.size() : cloud->size()); ++pid )
        {
            int index = indices?indices->indices[pid]:pid;
            labels_points[ (*labels)[index] ].push_back( index );
        }

    // display 300 lines
    float chance = 1000.f / lines.size();

    vector<vector<Scalar > > lines_ends;
    vector<int> calced_labels;
    if ( dim == 6 )
    {
        split( lines_ends, optimized_mask, lines, cloud, gap_threshold, &calced_labels );
        std::cout << "split finished" << std::endl;
    }

    std::vector<std::pair<int,PrimitiveT> > filtered; filtered.reserve( lines.size() );
    for ( size_t k = 0; k != lines.size(); ++k )
    {
        if ( !optimized_mask[k] && (rand()/static_cast<float>(RAND_MAX) > chance) ) continue;

        std::vector<MyPoint> minMax;
        int err = lines[k].getExtent( minMax
                                      , cloud
                                      , scale
                                      , labels ? &(labels_points[k])
                                               : indices ? &(indices->indices)
                                                         : NULL
                                      , true );
        if ( err != EXIT_SUCCESS )  { /*std::cout << "skipping "  << k << std::endl;*/ continue; }
        else
        {
            filtered.push_back( std::pair<int,PrimitiveT>(k,lines[k]) );
        }

        vector<pcl::PointXYZ> ps;
        if ( dim == 6 )
        {
            ps.push_back( am::util::pcl::asPointXYZ(minMax[0].getVector3fMap()
                          - .015f * (minMax[1].getVector3fMap() - minMax[0].getVector3fMap()) ) );
            ps.push_back( am::util::pcl::asPointXYZ(minMax[1].getVector3fMap()
                          - .05f * (minMax[0].getVector3fMap() - minMax[1].getVector3fMap()) ) );
        }
        else if ( dim == 4 )
        {
            pcl::PointXYZ tmp;
            for ( int corner_id = 0; corner_id != minMax.size(); ++corner_id )
            {
                tmp.x = minMax[corner_id].x; tmp.y = minMax[corner_id].y; tmp.z = minMax[corner_id].z;
                ps.push_back( tmp );
            }
        }
        else
            std::cerr << "[" << __func__ << "]: " << "It has to be lines or planes...\n";

        for ( int vpid = 0; vpid != 1; ++vpid )
        {
            PrimitiveT::draw( ps, vptrs[vpid], getLineName(k,vpid), 0., .3, 1. , vpids[vpid] );
        }

        if ( show_line_id )
            for ( int vpid = 0; vpid != 2; ++vpid )
                if ( vpid != 1 || optimized_mask[k] )
                    vptrs[vpid]->addText3D( am::util::sprintf("%d",k)
                                            , ps[0]
                            , 0.01
                            , 1., 0., 0.
                            , getLineTextName( k, vpid )
                            , vpids[vpid] );

        endpoints.insert( endpoints.end(), ps.begin(), ps.end() );
    }
    vptrs[0]->resetCamera(); vptrs[1]->resetCamera();
    vptrs[0]->spinOnce(); vptrs[1]->spinOnce();

    const int K = std::accumulate( optimized_mask.begin(), optimized_mask.end(), 0 );
    std::vector<Eigen::Vector3f> colours = am::util::nColoursEigen( K, 1.f, true );
    int colour_id = 0;

    // LINE RELATIONs
    std::ofstream f("line_angles.txt");
    //vector<float> ang_dists;
    vector<std::pair<float,float> > ang_dists;

    // show
    std::map<int,int> lineids_colourids;
    int               correct_count = 0;
    const int         vpid          = 1;
    std::vector<int>  finals;
    int               chosen_id     = 0;
    for ( size_t fid = 0; fid != filtered.size(); ++fid )
    {
        const int line_id = filtered[fid].first;

        // skip unchosen lines
        if ( !optimized_mask[line_id] )
        {
            ///bool suc = vptrs[vpid]->removeShape ( getLineName    (line_id,vpid), vpids[vpid] );
            //if ( !suc ) std::cout << "removed " << getLineName(line_id,vpid) << ((suc)?" OK":" NO") << std::endl;

            //vptrs[vpid]->removeText3D( getLineTextName(line_id,vpid), vpids[vpid] );
        }
        else // redraw chosen ones in colour
        {
            //bool ok = vptrs[vpid]->removeShape( getLineName(line_id,vpid), vpids[vpid] );
            //if ( !ok ) std::cerr << "could not find " << getLineName(line_id,vpid) << std::endl;

            vector<pcl::PointXYZ> ps;
            if ( dim == 6 )
            {
                //ps.insert( ps.end(), &(endpoints[fid*2]), &(endpoints[fid*2+1]) );
                if ( lines_ends[chosen_id].size() && lines_ends[chosen_id][0] < lines_ends[chosen_id][1] )
                {
                    for ( int i = 0; i != lines_ends[chosen_id].size(); i+= 2)
                    {
                        ps.push_back( pntToPCL<pcl::PointXYZ,Scalar>( lines[line_id].pos() + lines[line_id].dir() * lines_ends[chosen_id][0] ) );
                        ps.push_back( pntToPCL<pcl::PointXYZ,Scalar>( lines[line_id].pos() + lines[line_id].dir() * lines_ends[chosen_id][1] ) );

                        PrimitiveT::draw( ps
                                          , vptrs[vpid]
                                          , getLineName(line_id,vpid)
                                          , colours[colour_id][0], colours[colour_id][1], colours[colour_id][2]
                                          , vpids[vpid] );
                    }
                }
                else
                {
                    ps.push_back( endpoints[fid*2] );
                    ps.push_back( endpoints[fid*2+1] );

                    PrimitiveT::draw( ps
                                      , vptrs[vpid]
                                      , getLineName(line_id,vpid)
                                      , colours[colour_id][0], colours[colour_id][1], colours[colour_id][2]
                                      , vpids[vpid] );
                }
                ++chosen_id;

                if ( ps[1].getVector3fMap() != endpoints[fid*2+1].getVector3fMap() )
                     std::cerr << "wtf: "
                               << ps[1].x << ", " << ps[1].y << ", " << ps[1].z << "\t"
                               << endpoints[fid*2+1].x << ", " << endpoints[fid*2+1].y << ", " << endpoints[fid*2+1].z << std::endl;
//                vptrs[vpid]->addLine( endpoints[fid*2]
//                                      , endpoints[fid*2+1]
//                                      , colours[colour_id][0], colours[colour_id][1], colours[colour_id][2]
//                                      , getLineName(line_id,vpid)
//                                      , vpids[vpid] );

            }
            else if ( dim == 4 )
            {
                ps.insert( ps.end(), &(endpoints[fid*4]), &(endpoints[fid*4+4]) );
                PrimitiveT::draw( ps
                                  , vptrs[vpid]
                                  , getLineName(line_id,vpid)
                                  , colours[colour_id][0], colours[colour_id][1], colours[colour_id][2]
                                  , vpids[vpid] );
            }
            else
                std::cerr << "[" << __func__ << "]: " << "It has to be lines or planes...\n";

//            PrimitiveT::draw( ps
//                              , vptrs[vpid]
//                              , getLineName(line_id,vpid)
//                              , colours[colour_id][0], colours[colour_id][1], colours[colour_id][2]
//                              , vpids[vpid] );

            // draw weights
            Eigen::Vector3f centroid( Eigen::Vector3f::Zero() );
            for ( int c = 0; c != NCorners; ++c )
                centroid += endpoints[fid*NCorners+c].getVector3fMap();
            centroid /= static_cast<float>(NCorners);

            for ( size_t i = 0; i < finals.size(); ++i )
            {
                const int fid_i = finals[ i ];

                std::pair<float,float> score_angle = am::LocalAnalysis::matching( filtered[fid_i].second.dir()
                                                                                  , filtered[fid].second.dir()
                                                                                  , desired_angles );
                score_angle.second *= 180.f / M_PI;
                //float angle = score_angle.second;
                float angle = acos(filtered[fid_i].second.dir().dot( filtered[fid].second.dir() )) * 180.f / M_PI;
                // adjust nan-s
                if ( angle != angle ) angle = 0.f;

                // store for later
                ang_dists.push_back( std::pair<float,float>( score_angle.first, angle) );

                if ( 1e-5f > SimAnnOptProblem<PrimitivesT>::calcPairwiseTerm(filtered[fid_i].second,filtered[fid].second, desired_angles, 0.f) )
                {
                    continue;
                }

                // connect lines by centroids
                Eigen::Vector3f centroid_i( Eigen::Vector3f::Zero());
                for ( int c = 0; c != NCorners; ++c )
                    centroid_i += endpoints[fid_i*NCorners+c].getVector3fMap();
                centroid_i /= static_cast<float>(NCorners);

                vptrs[vpid]->addLine( smartgeometry::toPointXYZ(centroid), smartgeometry::toPointXYZ(centroid_i)
                                      , .45, .45, .45
                                      , "edge_" + getLineName(line_id,vpid) + getLineName(filtered[fid_i].first,vpid)
                                      , vpids[vpid] );

                // plot number
                float dev = std::min( fabs(90.f - angle), fabs(angle) ) / 45.f;
                vptrs[vpid]->addText3D( am::util::sprintf("%3.2f",angle)
                                        , smartgeometry::toPointXYZ((centroid_i+centroid)/2.f*(.9f + .1f * rand()/static_cast<float>(RAND_MAX)) )
                                        , 0.01
                                        , std::min(1.f,dev), std::max(1.f-dev,0.f), 0.
                                        , getLineTextName(line_id,vpid) + getLineTextName(filtered[fid_i].first,vpid)
                                        , vpids[vpid] );

            }

            lineids_colourids[ line_id ] = colour_id++;
            finals.push_back(fid);

            if ( real_mask && (*real_mask)[line_id] )
            {
                vptrs[vpid]->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3., getLineName(line_id,vpid), vpids[vpid] );
                    ++correct_count;
            }
        }
    } // show for

    // angles
    if ( show_angles )
    {
        std::cout<<"ang_dists:";for(size_t vi=0;vi!=ang_dists.size();++vi)std::cout<<ang_dists[vi].first<<" ";std::cout << "\n";
        std::sort( ang_dists.begin(), ang_dists.end(), [](std::pair<float,float> const& a, std::pair<float,float> const& b){ return a.first < b.first; } );
        for ( size_t vi=0; vi != ang_dists.size(); ++vi )
            f << ang_dists[vi].first * 180.f / M_PI << " "
              << ang_dists[vi].second << "\n";
        f.close();
        //int err = system( "gnuplot -e \"plot 'line_angles.txt' u 1 w linespoints, 'line_angles.txt' u 2 w points\" -p" );
        int err = system( "gnuplot -e \"plot 'line_angles.txt' u 1 w points\" -p" );
        if ( err != EXIT_SUCCESS ) std::cerr << "gnuplot returned " << err << std::endl;
    }

    // label colours
    pcl::PointCloud<MyPoint>::Ptr coloured_cloud( new pcl::PointCloud<MyPoint>() );
    pcl::copyPointCloud( *cloud, *coloured_cloud );
    if ( labels || calced_labels.size() )
    {
        auto const* p_labels = labels ? labels : &calced_labels;
        for ( size_t pid = 0; pid != coloured_cloud->size(); ++pid )
        {
            Eigen::Vector3f const& colour = colours[ lineids_colourids[(*p_labels)[pid]] ];
            am::util::pcl::setPointColor( coloured_cloud->at(pid),
                                          colour(0)*255.f, colour(1)*255.f, colour(2)*255.f );
        }

        std::string cloud_name = getCloudName(vpid);
        vptrs[vpid]->removePointCloud( cloud_name, vpids[vpid] );
        vptrs[vpid]->addPointCloud( coloured_cloud, cloud_name, vpids[vpid] );
        vptrs[vpid]->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3., cloud_name, vpids[vpid] );
    }

    // scale
    auto chosen = cloud->at(cloud->size()/2);
    vptrs[vpid]->addSphere( chosen, scale, "scale_sphere", vpid );
    vptrs[vpid]->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, 0.5, "scale_sphere", vpid );

    // complexity plot
    if ( 0 )
    {
//        int real_count = real_mask ? std::accumulate( real_mask->begin(), real_mask->end(), 0 ) : K;
        int opt_count = std::accumulate( optimized_mask.begin(), optimized_mask.end(), 0 );

        vptrs[vpid]->addText3D( am::util::sprintf("(opt) %d",opt_count/*,real_count*/) /*=?= %d (real)*/
                                , pcl::PointXYZ(0,0,0)
                                , 0.1
                                , 1.0, 0., 0.
                                , "success_rate"
                                , vpids[vpid] );
    }

    // show
    {
        vptrs[0]->spinOnce();
        if ( spin ) vptrs[vpid]->spin();
        else        vptrs[vpid]->spinOnce();
    }

    return EXIT_SUCCESS;
} // ... visualize

template <class PrimitivesT, typename Scalar> inline Scalar
am::GlobFit2::queryEnergy( Eigen::Matrix<Scalar,-1,1>                & energies
                           , MaskType                           const& mask
                           , PrimitivesT                        const& lines
                           , pcl::PointCloud<MyPoint>::ConstPtr        cloud
                           , Eigen::Matrix<Scalar,-1,1>         const& lambdas
                           , Scalar                             const  scale
                           , std::vector<Scalar>                const& desired_angles
                           , double                             const  threshold
                           , float                              const  trunc_at_pw_angle
                           , std::vector<int>                   const* indices
                           )
{
    SimAnnOptProblem<PrimitivesT> sa_problem( lines, cloud, threshold, indices, lambdas, scale, 65536, trunc_at_pw_angle );
    sa_problem.init( mask );
    sa_problem.calcEnergy( mask, desired_angles, &energies, true );
    return sa_problem.getEnergy( &desired_angles );
} // ... queryEnergy

#endif // __GF2_GLOBFIT2_HPP__
