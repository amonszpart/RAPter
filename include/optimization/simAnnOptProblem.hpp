#ifndef __GF2_SIMANNOPTPROBLEM_HPP__
#define __GF2_SIMANNOPTPROBLEM_HPP__

#include "optimization/simAnnOptProblem.h"
//#include "primitives/linePrimitive.h"

namespace am
{
    template <class PrimitivesT>
    SimAnnOptProblem<PrimitivesT>::SimAnnOptProblem( PrimitivesT             const& lines
                                        , pcl::PointCloud<MyPoint>::ConstPtr        cloud
                                        , double                             const  threshold
                                        , std::vector<int>                   const* indices
                                        , Eigen::VectorXf                    const& lambdas
                                        , float                                     scale
                                        , int                                       max_step_count
                                        , float                              const  trunc_pw_at_angle
                                        , PointPrimitiveDistanceCache<Scalar>     * point_prim_dist_cache )
        : _lines                ( lines             )
        , _configs_energies_capacity( 20            )
        , _lambdas              ( lambdas           )
        , _cloud                ( cloud             )
        , _threshold            ( threshold         )
        , _scale                ( scale             )
        , _energy_up_to_date    ( false             )
        , _step_count           ( 0                 )
        , _max_step_count       ( max_step_count    )
        , _smooth_scores        ( Eigen::Matrix<Scalar,-1,-1>::Zero(cloud->size(),cloud->size()) )
        , _subsample_lines_ratio( 1.f               )
        , _points_prims_dist_cache( point_prim_dist_cache )
        , _trunc_pw_at_angle    ( trunc_pw_at_angle )
        , _omp_thread_id        ( -1                )
        , _sensor               ( Eigen::Vector3f::Zero() )
        , _statistics           ( _max_step_count   )
    {
        if ( indices )
            _indices = *indices;

        // fill distances
        if ( lambdas(3) > 0.f )
            for ( size_t pid = 0; pid != _cloud->size()-1; ++pid )
                for ( size_t pid2 = pid+1; pid2 != _cloud->size(); ++pid2 )
                {
                    _smooth_scores(pid,pid2) = exp( -1.f * (_cloud->at(pid).getVector3fMap()-_cloud->at(pid2).getVector3fMap()).squaredNorm() / 25 );
                }

#if GF2_USE_CUDA
        // init gf2Cuda
        int dim = lines[0]().rows();
        std::vector<float> lines_data; lines_data.reserve( lines.size() * dim );
        for ( size_t i = 0; i != lines.size(); ++i ) for ( int d = 0; d != dim; ++d ) lines_data.push_back( lines[i]()(d) );
        _gf2_cuda.setLines( lines_data );
#endif
        std::cerr << "[" << __func__ << "]: " << "trunc_at: " << _trunc_pw_at_angle << " radians = " << _trunc_pw_at_angle * 180.f / M_PI << std::endl;
    }

    template <class PrimitivesT> typename SimAnnOptProblem<PrimitivesT>::Scalar
    SimAnnOptProblem<PrimitivesT>::calcDataTerms( pcl::PointCloud<MyPoint>::ConstPtr         cloud
                                                  , PrimitivesT                       const& lines
                                                  , MaskType                          const& line_mask
                                                  , std::vector<int>                  const* indices
                                                  , float                             const* scale
                                                  , std::vector<int>                       * labels
                                                  , PointPrimitiveDistanceCache<Scalar>    * cache
                                                  , Eigen::Matrix<float,3,1>          const* sensor )
    {
        const bool use_indices  = indices && indices->size();
        const int  N            = use_indices ? indices->size() : cloud->size();
        const int  L            = lines.size();

        // get line inliers
        std::vector<int>    points_lines ( N, -1      );
        std::vector<Scalar> min_distances_sqr( N, DBL_MAX );
        {
            // init sac model
            //pcl::SampleConsensusModelLine<MyPoint> sacline( cloud );
            //if ( indices && indices->size() )   sacline.setIndices( *indices );
            std::vector<Scalar> tmp_distances_sqr( N, DBL_MAX );

            // assignments - get closest line for each point
            for ( int line_id = 0; line_id != L; ++line_id )
            {
                // skip unchosen lines
                if ( !line_mask[line_id] ) continue;

                // cache hit?
                //if ( cache && cache->at(line_id).size() )   tmp_distances_sqr = cache->at( line_id );
                if ( cache && cache->contains(line_id) )   tmp_distances_sqr = cache->get( line_id );
                else
                {
                    // calc distances
                    //sacline.getDistancesToModel( lines[line_id](), tmp_distances );
                    if ( use_indices )
                        for ( size_t pid_id = 0; pid_id != N; ++pid_id )
                        {
                            tmp_distances_sqr[pid_id] = lines[line_id].point3Distance( cloud->at(indices->at(pid_id)).getVector3fMap() );
                            tmp_distances_sqr[pid_id] *= tmp_distances_sqr[pid_id];
                        }
                    else
                        for ( size_t pid    = 0; pid    != N; ++pid    )
                        {
                            tmp_distances_sqr[pid   ] = lines[line_id].point3Distance( cloud->at(pid                ).getVector3fMap() );
                            tmp_distances_sqr[pid   ] *= tmp_distances_sqr[pid];
                        }

                    // store in cache
                    //if ( cache )    cache->at(line_id) = tmp_distances_sqr;
                    if ( cache )    cache->push( line_id, tmp_distances_sqr );
                }

                // minimize distance for each point
                for ( int pnt_id = 0; pnt_id != N; ++pnt_id )
                    if ( tmp_distances_sqr[pnt_id] < min_distances_sqr[pnt_id] )
                    {
                        min_distances_sqr[pnt_id] = tmp_distances_sqr[pnt_id];
                        points_lines     [pnt_id] = line_id;
                    }
            } // for line_id
        } // get line inliers

        if ( labels )
            *labels = points_lines;

        // return squared distance
        //return std::inner_product( min_distances.begin(), min_distances.end(), min_distances.begin(), 0.f );
//        struct accumFunctor
//        {
//                float _scaleSqr;
//                accumFunctor ( float scale ) : _scaleSqr( scale * scale ) {}
//                float operator()( float i, float const& x )
//                {
//                    float dist =
//                    return i + x * x / _scaleSqr;
//                }
//        };
        //return std::accumulate( min_distances.begin(), min_distances.end(), 0.f, accumFunctor(*scale) );

        const Scalar scaleSqr = (*scale) * (*scale);
        Scalar ret( 0 );
        for ( size_t pid = 0; pid != min_distances_sqr.size(); ++pid )
        {
            Scalar distSqr = min_distances_sqr[pid];
//            if ( sensor )
//            {
//                Scalar sensor_dist = (cloud->at(pid).getVector3fMap() - *sensor).norm();
//                distSqr /= sensor_dist * sensor_dist;
//            }
            ret += distSqr;
        }
        ret /= scaleSqr * static_cast<Scalar>( cloud->size() );

        return ret;
    }

    template <class PrimitivesT> inline typename SimAnnOptProblem<PrimitivesT>::Scalar
    SimAnnOptProblem<PrimitivesT>::calcPairwiseTerm( PrimitiveT                                      const &l1
                                                    , PrimitiveT                                     const &l2
                                                    , std::vector<Scalar>                            const &desired_angles
                                                    , Scalar const trunc_at
                                                    , bool                                           const verbose
                                                   )
    {

        Scalar angle = acos( l1.dir().dot(l2.dir()) );
        // nan
        if ( angle != angle )   angle = 0.f;
        if ( verbose ) std::cout << "angle " << angle * 180.f / M_PI;

        Scalar min_ang_diff = FLT_MAX, tmp;
        for ( int ang_id = 0; ang_id != desired_angles.size(); ++ang_id )
            if ( (tmp=fabs(angle - desired_angles[ang_id])) < min_ang_diff )
                min_ang_diff = tmp;

        const Scalar trunc_at2 = trunc_at * 2.f;
        Scalar ret = FLT_MAX;
        if ( min_ang_diff < trunc_at )
        {
            //tmp = min_ang_diff / trunc_at;
//            ret = tmp * tmp * tmp;// * tmp * tmp * tmp
//            ret *= ret;
            if ( min_ang_diff < 1e-6f )
                ret = 0.f;
            else
                ret = std::exp( 0.1f + 50.f * min_ang_diff ) - static_cast<Scalar>(1);
        }
//        else if ( min_ang_diff < trunc_at2)
//        {
//            tmp =  (min_ang_diff - trunc_at2) / trunc_at;
//            ret = tmp * tmp * tmp;// * tmp * tmp * tmp;
//            ret *= ret;
//        }
        else
            ret = 0.f;

        if ( !std::isfinite(ret) )
        {
            std::cerr << "ret: " << ret << ", for lines " << l1().transpose()
                                      << ", " << l2().transpose()
                                      << ", angle: " << angle
                                      << ", min_ang_diff: " << min_ang_diff
                                      << ", exp(min_ang_diff): " << exp(min_ang_diff)
                                      << ", exp( 0.1f + min_ang_diff): " << exp(0.1f + min_ang_diff )
                                      << ", exp( 0.1f + 500.f * min_ang_diff): " << exp(0.1f + 500.f * min_ang_diff )
                                      << std::endl;

            ret = FLT_MAX;
        }
        if ( verbose ) std::cout << ", mined angle " << min_ang_diff * 180.f / M_PI << " <? " << trunc_at * 180.f / M_PI << " = " << ret << std::endl;

        return ret;
    }

    template <class PrimitivesT> typename SimAnnOptProblem<PrimitivesT>::Scalar
    SimAnnOptProblem<PrimitivesT>::calcEnergy( MaskType                     const& line_mask
                                               , std::vector<Scalar>        const& desired_angles
                                               , Eigen::Matrix<Scalar,-1,1>      * energies
                                               , bool                       const  verbose ) const
    {
        typedef typename PrimitivesT::value_type PrimitiveT;

        //const int N = _cloud->size();
        const int K = _lines.size();
        std::vector<int> labels;

        Eigen::Matrix<Scalar,-1,1> curr_energies( Eigen::Matrix<Scalar,4,1>::Zero() );

        // E_complexity
        {
            curr_energies( 1 ) = std::accumulate( line_mask.begin(), line_mask.end(), 0 );
            //std::cout << "complexity term is " << e_complexity[comb_id] << std::endl;
            //if ( curr_energies(1) <= 0.f )   return curr_energies(0);
        }

        // E_data
        //if ( _lambdas(0) > 0.f )
        {
            curr_energies(0) = calcDataTerms( /*        cloud: */ _cloud
                                              , /*      lines: */ _lines
                                              , /*  curr_mask: */ line_mask
                                              , /*  p_indices: */ ((_indices.size() > 0) ? &_indices : NULL)
                                              , /*      scale: */ &_scale
                                              , /* out_labels: */ &labels
                                              , /*      cache: */ _points_prims_dist_cache
                                              , &_sensor );
//            std::cout << "curr_energies(0): " << curr_energies(0) << std::endl;
        } // E_data

        // E_pw
        //if ( _lambdas(2) > 0.f )
        {
#if GF2_USE_CUDA
            //float pwEnergyGPU = 0.f;
            _gf2_cuda.getEnergyOfHostMaskDense( curr_energies(2), line_mask );
#else
            int count = 0;
            // for every i-j pair
            for ( int i = 0; i != K-1; ++i )
            {
                // skip unchosen lines
                if ( !line_mask[i] ) continue;
                // cache
                PrimitiveT const& line0 = _lines[ i ];

                for ( int j = i+1; j != K; ++j )
                {
                    if ( !line_mask[j] ) continue;

                    if ( verbose ) std::cout << "comparing lines " << i << "<->" << j << " ";

                    curr_energies(2) += calcPairwiseTerm( line0, _lines[j], desired_angles, _trunc_pw_at_angle /*0.4f*/, verbose ); // trunc at 22.5 degrees
                    ++count;
                }
            }
//            std::cout << "curr_energies(2): " << curr_energies(2);
            if ( count )
                curr_energies(2) /= static_cast<Scalar>( count );
//            std::cout << " --> " << curr_energies(2) << std::endl;
            //if ( pwEnergyGPU == curr_energies(2) ) std::cout << "curr_energies(2) " << curr_energies(2) << " == " << pwEnergyGPU << " energyGPU\n";
            //else                         std::cerr << "curr_energies(2) " << curr_energies(2) << " == " << pwEnergyGPU << " energyGPU\n";
#endif
        } // E_pw

        // E_smooth labeling
//        if ( _lambdas(3) > 0.f )
//        {
//            curr_energies(3) = 0.f;
//            for ( size_t pid = 0; pid != _cloud->size()-1; ++pid )
//                for ( size_t pid2 = pid+1; pid2 != _cloud->size(); ++pid2 )
//                {
//                    if ( labels[pid] != labels[pid2] )
//                        //curr_energies(3) += exp( -1.f * (_cloud->at(pid).getVector3fMap()-_cloud->at(pid2).getVector3fMap()).squaredNorm() / 25 );
//                        curr_energies(3) += _smooth_scores( pid, pid2 );
//                }
//        }

        // final E
        if ( energies )     *energies = curr_energies;

        return _lambdas.dot( curr_energies );
    }

    template <class PrimitivesT> typename SimAnnOptProblem<PrimitivesT>::Scalar
    SimAnnOptProblem<PrimitivesT>::getEnergy( void const* desired_angles )
    {
        if ( !desired_angles ) { std::cerr << "[" << __func__ << "]: " << "need desired angles here..." << std::endl; return FLT_MAX; }
        // return if cached
        if ( _energy_up_to_date )   return _configs_energies.front().second;

        // selection
        MaskType const& line_mask = _configs_energies.front().first;

        // final E

        _configs_energies.front().second = calcEnergy( line_mask, *static_cast<std::vector<Scalar>const*>(desired_angles) );
        _energy_up_to_date               = true;

        return _configs_energies.front().second;
    }

    template <class PrimitivesT> int
    SimAnnOptProblem<PrimitivesT>::init( std::vector<int> const& config )
    {
        _configs_energies.push_front( std::pair<MaskType,Scalar>(config,-1) );

        return EXIT_SUCCESS;
    }

    template <class PrimitivesT> inline std::vector<int>
    sampleCluster( int                          const  sample_count
                   , std::vector<int>           const& cluster_lines
                   , PrimitivesT                const& lines
                   , int                        const  max_retries             //= 20
                   , int                        const* representative_id
                )
    {
        const static int RAND_MAX_2 = RAND_MAX/2;
        // clusters nowadays don't have multiple directions
        if ( sample_count != 1 )
        {
            std::cerr << "[" << __func__ << "]: " << "sample_count != 1...exiting\\n";
            return std::vector<int>();
        }
        std::vector<int> lids(1);

        // return representative direction
        if ( representative_id && rand() < RAND_MAX_2 )
            if ( (static_cast<size_t>(*representative_id) < lines.size()) && (*representative_id >= 0) )
            {
                lids[0] = *representative_id;
                return lids;
            }
            else
                std::cerr << "[" << __func__ << "]: " << "representative_id provided, but invalid..." << *representative_id << std::endl;
        else
            lids[0] = cluster_lines[ rand() % cluster_lines.size() ];
#if 0
        // do until perpendicular
        bool retry       = false;
        int  retry_count = 0;
        do
        {
            // reset
            retry = false;
            // sample
            for ( int i = 0; i != sample_count; ++i )   lids[i] = cluster_lines[ rand() % cluster_lines.size() ];

            // check perpendicularity
//            if ( !retry && (sample_count == 2) )
//            {
//                float cos_angle = fabs( lines[lids[0]].dir().dot( lines[lids[1]].dir()) );
//                retry |= (cos_angle > 0.5f);
//            }
        } while ( retry && ((max_retries < 1) || (retry_count++ < max_retries)) );

        // couldn't find it, skip one
        if ( retry_count >= max_retries )
            lids.resize(1);
#endif

        return lids;
    }

    /**
     * @brief clusteredInit Wiser initialisation of the optimisation using clustering information.
     * @param clustering    Information about line memberships in clusters.
     * @param config        Contains a selection mask of clusters to select cluestering.directions_in_clusters # of lines from.
     * @return
     */
    template <class PrimitivesT> int
    SimAnnOptProblem<PrimitivesT>::initWithClusters( PrimitiveClustering<PrimitivesT> const& clustering
                                                     , MaskType                       const& config
                                                     , int                            const  K         )
    {
        if ( clustering.directions_in_cluster > 2 )
        {
            std::cerr << "[" << __func__ << "]: " << "clustering.directions_in_cluster > 2...exiting\n";
            return EXIT_FAILURE;
        }

        //int              k                         = 0; // how many lines selected up til now
        MaskType         mask( _lines.size(), 0 );      // line selection mask (dense)
        std::set   <int> line_ids;                      // mask mirror
        std::vector<int> lids;                          // tmp storage for lines sampled from a cluster

        int full_retries = 0;
        do
        {
            // count how many we have
            //k = std::accumulate( mask.begin(), mask.end(), 0 );

            // for each chosen cluster coded in config
            for ( size_t cid = 0; (cid != config.size()) && (line_ids.size() < K); ++cid )
            {
                // skip, if not chosen
                if ( !full_retries && !config[cid] ) continue;

                // debug
#               pragma omp critical
                {
                    if ( !config[cid] )
                        if ( !_omp_thread_id ) std::cout << "full_retries: " << full_retries << ", config[i]: " << config[cid] << std::endl;
                }

                int     lids_unique = 0; // count unique new samples
                size_t  tmp_retries = 0; // cerr, if stuck
                do
                {
                    // reset
                    lids_unique = 0;

                    // sample
                    lids = sampleCluster( clustering.directions_in_cluster
                                          , clustering.clusters_lines[ config[cid] ]
                                          , _lines
                                          , 20
                                          , /* representative_id: */ NULL ); // we want different lines on every init
                    // count unique
                    for ( size_t lid_id = 0; lid_id != lids.size(); ++lid_id )
                        lids_unique += (line_ids.find( lids[lid_id] ) == line_ids.end());

#                   pragma omp critical
                    {
                        // STUCK check
                        if ( ++tmp_retries > clustering.clusters_lines[ config[cid] ].size() )
                        {
                            std::cerr << "[" << __func__ << "][" << _omp_thread_id << "]: " << "STUCK!!"
                                      << "size: " << static_cast<int>(clustering.clusters_lines[ config[cid] ].size()) << "/" << K
                                                                                                                     << std::endl;

                            std::cerr<<"clustering.clusters_lines[config[i]=="<<config[cid]<<"]:";for(size_t vi=0;vi!=clustering.clusters_lines[config[cid]].size();++vi)std::cout<<clustering.clusters_lines[config[cid]][vi]<<" ";std::cout << "\\n";

                            std::cerr << "line_ids:";
                            for ( auto it = line_ids.begin(); it != line_ids.end(); ++it )
                                std::cerr<< *it <<" ";
                            std::cerr << "\n";

                            std::cerr << "lids: " << lids[0];
                            if (lids.size() > 1) std::cerr << lids[1];
                            std::cerr << std::endl;

                            //if ( tmp_retries > static_cast<int>(clustering.clusters_lines[ config[i] ].size()) * 2 )
                            //    return EXIT_FAILURE;
                        }
                    }

                } while (    (lids_unique != clustering.directions_in_cluster)
                          && (tmp_retries <  clustering.clusters_lines[config[cid]].size()) );

                for ( size_t lid_id = 0; (lid_id != lids.size()) && (line_ids.size() < K); ++lid_id )
                {
                    //if ( !mask[ lids[lid_id] ] ) ++k;

                    line_ids.insert( lids[lid_id] );
                    mask[ lids[lid_id] ] = 1;
                } // ... for lids
            } // ... for configs

            if ( line_ids.size() < config.size() )
                while ( line_ids.size() < K )
                {
                    int line_id = rand() % mask.size();
                    if ( (line_ids.find( line_id ) == line_ids.end()) && !mask[line_id] )
                    {
                        line_ids.insert( line_id );
                        mask[line_id] = 1;
                    }
                }
        } while ( (line_ids.size() < K) && (full_retries++ < 10) );

        if ( K != std::accumulate( mask.begin(), mask.end(), 0 ) )
        {
            std::cerr << "[" << __func__ << "]" << "[" << _omp_thread_id << "]: "
                      << "maskK: " << std::accumulate( mask.begin(), mask.end(), 0 )
                      << "/" << K << std::endl;
            return EXIT_FAILURE;
        }

        this->init( mask );

        _clustering = clustering;
        return EXIT_SUCCESS;
    }

    template <class PrimitivesT> bool
    SimAnnOptProblem<PrimitivesT>::hasNext()
    {
        return ( _step_count + 1 < _max_step_count );
    }

    // 1. select best candidate to remove
    template <class PrimitivesT, typename Scalar = typename SimAnnOptProblem<PrimitivesT>::Scalar> inline std::pair<int,Scalar>
    selectLineToRemove( MaskType                     const& mask
                        , PrimitivesT                const& lines
                        , std::vector<Scalar>        const& desired_angles
                        , float                      const  one_over_current_step_count
                        , SimAnnOptProblem<PrimitivesT> const& sa    )
    {
        using std::vector;
        using std::pair;

        pair<int,Scalar> min_e( -1, FLT_MAX );

        MaskType tmp_mask = mask;
        for ( size_t lid = 0; lid != lines.size(); ++lid )
        {
            // look for 1-s
            if ( !mask[lid] ) continue;

            // stochastic speedup
            //if ( (_subsample_lines_ratio < 1.f) && (rand()/static_cast<float>(RAND_MAX) > _subsample_lines_ratio) ) continue;

            // flip
            tmp_mask[ lid ] = 0;

            // minimize energy
            Scalar e = sa.calcEnergy( tmp_mask, desired_angles );
            if ( !std::isfinite(e) )                e = FLT_MAX;
            // Simulated Annealing: keep high energy answer with 1/stepid probability
            //std::cout << "e: " << e << ", min_e.second: " << min_e.second << ", min_e.first: " << min_e.first << std::endl;
            if ( (e < min_e.second) && ((min_e.first < 0) || (rand()/static_cast<float>(RAND_MAX) > one_over_current_step_count)) )
            {
                min_e.second = e;
                min_e.first  = lid;
            }

            // restore
            tmp_mask[ lid ] = 1;
        } // ... for lid

        //std::cout << "[" << __func__ << "]: " << "returning " << min_e.first << std::endl;
        return min_e;
    }

    // 2. select cluster, to select best candidate from
    template <class PrimitivesT, class Scalar = typename SimAnnOptProblem<PrimitivesT>::Scalar> inline std::pair<int,Scalar>
    selectClusterToAddFrom( MaskType                         const& mask
                            , int                            const  line_id_to_remove
                            , std::vector<std::vector<int> > const& clusters_lines
                            , std::vector<int>               const& lines_clusters
                            , int                            const  directions_in_cluster
                            , int                            const  first_c
                            , std::vector<int>               const& cluster_representatives
                            , float                          const  one_over_current_step_count
                            , std::vector<int>               const& cluster_id_tabus
                            , PrimitivesT                    const& lines
                            , std::vector<Scalar>            const& desired_angles
                            , SimAnnOptProblem<PrimitivesT>  const& sa    )
    {
        using std::pair;
        using std::vector;
        //typedef typename SimAnnOptProblem<PrimitivesT>::Scalar Scalar;

        pair<int,Scalar> min_e( -1, FLT_MAX ); // <cid,min_e>
        // copy current
        MaskType tmp_mask = mask;

        // exclude worst candidate from step 1.
        tmp_mask[ line_id_to_remove ] = 0;

        // apply repulsion rules: don't allow this cluster, if an l_i is already chosen
//        std::vector<int> forbidden_clusters;
//        {
//            // over all lines
//            for ( size_t lid = 0; lid != lines.size(); ++lid ) // iterate over chosen lines
//            {
//                // chosen lines needed
//                if ( !tmp_mask[lid] ) continue;

//                // lids < first_c repulses it's cluster
//                if ( static_cast<int>(lid) < first_c ) forbidden_clusters.push_back( lines_clusters[ lid ] );
//            }
//        }

        // include one new candidate
        for ( size_t cid = 0; cid != clusters_lines.size(); ++cid )
        {
            // skip failed clusters
            if ( std::find(cluster_id_tabus.begin(), cluster_id_tabus.end(), cid) != cluster_id_tabus.end() ) continue;

            // don't allow forbidden clusters with temperature probability (they still contain l_i lines, so should be allowed)
//            if (    (std::find( forbidden_clusters.begin(), forbidden_clusters.end(), cid ) != forbidden_clusters.end())
//                 && (rand()/static_cast<float>(RAND_MAX) > one_over_current_step_count)                                  )
//                continue;


            // we might need to resample lines from the cluster to find at least one good candidate
            int retry_count = 0;
            bool resample = false;
            do
            {
                // reset
                resample = false;

                // get two lines from cluster
                if ( cluster_representatives.size() <= cid ) std::cerr << "[" << __func__ << "]: " << "WTF...cluster_representatives.size() <= cid: " << cluster_representatives.size() << " <= " << cid << std::endl;

                vector<int> lids = sampleCluster( /*   how many: */ directions_in_cluster
                                                  , /* line_ids: */ clusters_lines[cid]
                                                  , /*    lines: */ lines
                                                  , 20
                                                  , &(cluster_representatives[cid]) );

                // this usually goes to 2 for clusters of perpendicular lines
                for ( size_t lid_id = 0; lid_id != lids.size(); ++lid_id )
                {
                    // look for 0-s to flip to 1-s
                    if ( tmp_mask[ lids[lid_id] ] ) { resample = true ; continue; }
                    else                            { resample = false;           }

                    // flip ON
                    tmp_mask[ lids[lid_id] ] = 1;

                    // minimize energy
                    Scalar e = sa.calcEnergy( tmp_mask, desired_angles );
                    //std::cout << "1/s: " << one_over_current_step_count << ": " << std::endl;
                    if ( (min_e.first < 0) || ((e < min_e.second) && ( (rand()/static_cast<float>(RAND_MAX)) > one_over_current_step_count)) )
                    {
                        min_e.second = e;
                        min_e.first  = cid;
                    }

                    // flip OFF to restore for next cycle
                    tmp_mask[ lids[lid_id] ] = 0;
                } // ... for all samples of cluster
            } while ( resample && (retry_count++ < 20) ); // until at least one calcEnergy was run
        } // ... for each cluster

        return min_e;
    }

    // 3. select line to add from cluster selected
    template <class PrimitivesT, class Scalar = typename SimAnnOptProblem<PrimitivesT>::Scalar> inline std::pair<int,Scalar>
    selectLineToAddFromCluster( MaskType                         const& mask
                                , int                            const  line_id_to_remove
                                , int                            const  cluster_id
                                , std::vector<std::vector<int> > const& clusters_lines
                                , std::vector<int>               const& lines_clusters
                                , int                            const  first_c
                                , float                          const  scale
                                , PrimitivesT                    const& lines
                                , std::vector<Scalar>            const& desired_angles
                                , SimAnnOptProblem<PrimitivesT>  const& sa
                                )
    {
        //typedef typename SimAnnOptProblem<PrimitivesT>::Scalar Scalar;
        using std::pair;
        using std::vector;

        const int max_retries = 1; // (_subsample_lines_ratio < 1.f) ? (1.f / _subsample_lines_ratio) : 1;

        pair<int,Scalar> min_e2( -1, FLT_MAX );
        // copy current
        MaskType tmp_mask = mask;
        // exclude worst candidate from step 1.
        tmp_mask[ line_id_to_remove ] = 0;

        //const Eigen::Vector3f x0( (Eigen::Vector3f()<<0.f,0.f,0.f).finished() );

        // apply repulsion rules:
        std::vector<int> clusters_of_local_lines, clusters_of_consensus_lines;
        {
            // over all lines
            for ( size_t lid = 0; lid != lines.size(); ++lid ) // iterate over chosen lines
            {
                // chosen lines needed
                if ( !tmp_mask[lid] ) continue;

                // lids < first_c repulses it's cluster
                if ( static_cast<int>(lid) < first_c ) clusters_of_local_lines    .push_back( lines_clusters[ lid ] );
                else                                   clusters_of_consensus_lines.push_back( lines_clusters[ lid ] );
            }
        }

        // include one new candidate
        int retries = 0;
        while ( (min_e2.first < 0) && (retries++ < max_retries) )
        {
            // for all lines in cluster
            for ( size_t lid_id = 0; lid_id != clusters_lines[cluster_id].size(); ++lid_id )
            {
                const int lid = clusters_lines[cluster_id][ lid_id ];
                const int cid = lines_clusters[lid];

                // look for 0-s and skip candidate from step 1.
                if ( mask[lid] || (lid == line_id_to_remove) ) continue;

                // check, if forbidden connection
                if ( lid < first_c) // local approx line is not compatible with it's own cluster's consensus directions
                {   // local line
                    if ( std::find(clusters_of_consensus_lines.begin(), clusters_of_consensus_lines.end(), cid) != clusters_of_consensus_lines.end() )
                        continue;
                }
                else
                {   // consensus line
                    if ( std::find(clusters_of_local_lines.begin(), clusters_of_local_lines.end(), cid) != clusters_of_local_lines.end() )
                        continue;
                }
#if 0
                // skip lid_id, if too close in d space to any line in mask
                {
                    // d of current candidate
                    float r0 = (x0-lines[lid]().template segment<3>(0)).cross(
                                   lines[lid]().template segment<3>(3)).norm();

                    // skip, if too close
                    bool skip = false;
                    // for all lines selected by mask
                    for ( int i = 0; i != static_cast<int>(lines.size()) && !skip; ++i )
                    {
                        // only selected lines
                        if ( !tmp_mask[i] ) continue;

                        // d of line in mask
                        float r = (x0-lines[i]().template segment<3>(0)).cross(
                                      lines[i]().template segment<3>(3)).norm();

                        // check, if too close
                        skip |= ( fabs(r - r0) < scale );

                    }
                    if ( skip ) continue;
                }
#endif

                // stochastic speedup
                //if ( (_subsample_lines_ratio < 1.f) && (rand()/static_cast<float>(RAND_MAX) > _subsample_lines_ratio) ) continue;

                // flip
                tmp_mask[ lid ] = 1;

                // minimize energy
                Scalar e = sa.calcEnergy( tmp_mask, desired_angles );
                if ( e < min_e2.second )
                {
                    min_e2.second = e;
                    min_e2.first = lid;
                }

                // restore
                tmp_mask[ lid ] = 0;

            } // ... for lid
        } // ... while not chosen at least one

        return min_e2;
    }

    /**
     * @brief SimAnnOptProblem::next
     * doublegreedy: 1: select the worst line from the current config, remove it. 2. select the best new candidate from all the other lines, add it.
     * @return
     */
    template <class PrimitivesT> int
    SimAnnOptProblem<PrimitivesT>::next( void const* p_desired_angles )
    {
        using std::vector;
        using std::pair;

        // sanity check
        if ( !_configs_energies.size() )
        {
            std::cerr << "[" << __func__ << "] " << "please init before calling..." << std::endl;
            return EXIT_FAILURE;
        }

        if ( !p_desired_angles )
        {
            std::cerr << "[" << __func__ << "]: " << "!p_desired_angles...exiting\n";
            return EXIT_FAILURE;
        }

        // invalidate status
        _energy_up_to_date = false;

        // update
        MaskType mask = _configs_energies.front().first;

        // reverse jump
        {
            const unsigned    backLimit   = _configs_energies_capacity / 2;
            if ( (_configs_energies.size() > backLimit) && (_step_count > _max_step_count / 10.f) ) // only do it, if already enough ahead
            {
                const Scalar curr_energy = _configs_energies.front().second;
                for ( typename HistoryT::iterator it = _configs_energies.begin() + backLimit; (it != _configs_energies.end()); ++it )
                {
                    float prob = (curr_energy-it->second) / (it->second+curr_energy) * 2.f;

                    if ( prob < 0.f )   continue;

                    if ( (rand()/static_cast<float>(RAND_MAX)) < prob )
                    {
                        mask = it->first;
                        _statistics.addReverseJump( std::pair<int,Scalar>(_step_count, prob) );
                        _configs_energies.erase( _configs_energies.begin(), it );
                        break;
                    }
                }
            }
        }

        const int   K                   = std::accumulate( mask.begin(), mask.end(), 0 );
        const float one_over_step_count = std::min(0.5f,1.f / (_step_count+1.f));

        pair<int,Scalar> min_line_id_2rem   (  0, FLT_MAX );
        pair<int,Scalar> min_cluster_id_2add( -1, FLT_MAX ); // <cid,min_e>
        pair<int,Scalar> min_line_id_2add   ( -1, FLT_MAX ); // <lid,min_e>
        int retry_count = 0; const int max_retry_count = 10; //std::accumulate(mask.begin(),mask.end(),0);
        do
        {
            // 1. select best candidate to remove
            min_line_id_2rem = selectLineToRemove( mask, _lines, *static_cast<std::vector<Scalar> const*>(p_desired_angles), one_over_step_count, *this );

            // failsafes
            {
                if ( min_line_id_2rem.first < 0 )
                {
                    std::cerr << "[" << _omp_thread_id << "] " << "first not ok, retrying..." << std::endl;
                    ++retry_count;
                    continue;
                }
                if ( std::accumulate( mask.begin(), mask.end(), 0 ) != K )
                {
                    std::cerr << "[" << _omp_thread_id << "] " << "K not ok after first step in it " << _step_count << std::endl;
                }
                if ( !mask[min_line_id_2rem.first] ) std::cerr << "ooooooo" << std::endl;
            }

            vector<int> tabu_cluster_ids;
            do
            {
                // 2. select cluster, to select best candidate from
                //std::cout << "1/s: " << one_over_step_count << ": " << _step_count << std::endl;
                min_cluster_id_2add = selectClusterToAddFrom( mask
                                                              , min_line_id_2rem.first
                                                              , _clustering.clusters_lines
                                                              , _clustering.lines_clusters
                                                              , _clustering.directions_in_cluster
                                                              , _clustering.first_c
                                                              , _clustering.cluster_representatives
                                                              , one_over_step_count
                                                              , tabu_cluster_ids
                                                              , _lines
                                                              , *static_cast<std::vector<Scalar> const*>(p_desired_angles)
                                                              , *this );
                // safety break
                if ( min_cluster_id_2add.first < 0 )
                {
                    std::cerr << "[" << __func__ << "]: " << "min_e_cluster_add.first < 0...exiting\n";
                    return EXIT_FAILURE;
                }

                // 3. select best candidate to add instead
                min_line_id_2add = selectLineToAddFromCluster( mask
                                                               , min_line_id_2rem.first
                                                               , min_cluster_id_2add.first
                                                               , _clustering.clusters_lines
                                                               , _clustering.lines_clusters
                                                               , _clustering.first_c
                                                               , _scale
                                                               , _lines
                                                               , *static_cast<std::vector<Scalar> const*>(p_desired_angles)
                                                               , *this );

                // prevent this from happening
                if ( min_line_id_2add.first == -1 )
                {
                    tabu_cluster_ids.push_back( min_cluster_id_2add.first );
                }
            }
            while ( (min_line_id_2add.first == -1) && (tabu_cluster_ids.size() < _clustering.clusters_lines.size()) );

            // failsafes
            {
                if ( std::accumulate( mask.begin(), mask.end(), 0 ) != K )
                {
                    std::cerr << "[" << _omp_thread_id << "] " << "K not ok after all steps in it " << _step_count << std::endl;
                }
                if ( (min_line_id_2rem.first < 0)
                     || (min_cluster_id_2add.first < 0)
                     || (min_line_id_2add.first < 0) )
                {
                    std::cerr << "[" << _omp_thread_id << "] " << "Not ok..."
                              << min_line_id_2rem.first << ", "
                              << min_cluster_id_2add.first << ", "
                              << min_line_id_2add.first
                              << std::endl;
                    ++retry_count;
                    continue;
                }
            }
        }
        while ( (min_line_id_2add.first == -1) && (retry_count++ < max_retry_count) );

        if ( (min_line_id_2add.first == -1) )
        {
            std::cerr << " early termination..." << std::endl;
            _step_count = _max_step_count-1;
            return EXIT_FAILURE;
        }

        // apply changes
        {
            if ( !mask[ min_line_id_2rem.first ] ) std::cerr << "wtf..min_line_id_2rem not 1..." << std::endl;
            mask[ min_line_id_2rem.first ] = 0;
            if ( mask[ min_line_id_2add.first ] ) std::cerr << "wtf..min_line_id_2rem not 0..." << std::endl;
            mask[ min_line_id_2add.first ] = 1;
        }

        {
            if ( std::accumulate( mask.begin(), mask.end(), 0 ) != K )
                std::cerr << "[" << __func__ << "]["<<_omp_thread_id<<"]: " << "step 3 allowed " << std::accumulate( mask.begin(), mask.end(), 0 ) << " instead of " << K << std::endl;
        }
#if 1
        // 4. break loops
        {
            //int desiredK = std::accumulate( mask.begin(), mask.end(), 0 );
            int hammingSum = 0;
            // for every previous config in history
            //std::cout << "hamming: ";
            for ( auto it = _configs_energies.begin()+1; it != _configs_energies.end(); ++it )
            {
                // calculate hamming distance
                int dist = 0;
                for ( size_t lid = 0; lid != it->first.size(); ++lid )
                    dist += (mask[lid] != (it->first)[lid]);
                //std::cout << dist << " ";
                hammingSum += dist;
            }
            //float prob = std::max(0.f, 1.f - hammingSum / ( (float)(_configs_energies.size() - 1) * _configs_energies.size()) );
            float prob = std::max(0.f, 1.f - hammingSum / ((_configs_energies.size() - 1) * 2.f) );
            if ( rand()/static_cast<float>(RAND_MAX) < prob )
            {
                // look for a 1 and a 0
                int nids[2] = {0,0};
                do
                {
                    nids[0] = rand() % _lines.size();
                    nids[1] = rand() % _lines.size();
                } while ( (mask[nids[0]] + mask[nids[1]] != 1) || (nids[0] == nids[1]) );

                assert( nids[0] != nids[1] );
                mask[nids[0]] = !mask[nids[0]];
                mask[nids[1]] = !mask[nids[1]];

                min_line_id_2add.second = calcEnergy( mask, *static_cast<std::vector<Scalar> const*>(p_desired_angles) );

                _statistics.addHammingHit( std::pair<int,Scalar>(_step_count, prob) );
            }
        }
#endif

        // store result
        _configs_energies.push_front( std::pair<MaskType,Scalar>(mask, min_line_id_2add.second) );
        _energy_up_to_date = true;

        // limit deque
        if ( _configs_energies.size() > _configs_energies_capacity )    _configs_energies.pop_back();

        ++_step_count;

        // debug
        if ( !(_step_count % 100) )
//        if ( _step_count == _max_step_count-1 )
            std::cout << "[" << __func__ << "]" << "[" << _omp_thread_id << "]: "
                      << _step_count << "/" << _max_step_count << std::endl;

        return EXIT_SUCCESS;
    }

    template <class PrimitivesT> inline int
    SimAnnOptProblem<PrimitivesT>::stepBack()
    {
        _configs_energies.pop_front();
        _energy_up_to_date = true;

        return EXIT_SUCCESS;
    }

    template <class PrimitivesT> inline MaskType
    SimAnnOptProblem<PrimitivesT>::getConfig()
    {
        return _configs_energies.front().first;
    }

    template <class PrimitivesT> inline SimAnnStats
    SimAnnOptProblem<PrimitivesT>::report() const
    {
        std::cout << "statistics[" << _omp_thread_id << "]:\n"
                  << "\thamminghits: " << _statistics.hamming_hits.size() << "/" << _step_count+1 << " = " << _statistics.hamming_hits.size() / static_cast<float>( _step_count+1 ) * 100.f << " %\n"
                  << "\trevJumps: " << _statistics.rev_jumps.size() << "/" << _step_count+1 << " = " << _statistics.rev_jumps.size() / static_cast<float>( _step_count+1 ) * 100.f << " %\n";

        return _statistics;
    }

    template <class PrimitivesT> inline std::vector<MaskType>
    SimAnnOptProblem<PrimitivesT>::powerSet( int K )
    {
        if ( K == 1 )
            return {{0},{1}};

        std::vector< std::vector<int> > prev = powerSet( K-1 );

        std::vector< std::vector<int> > out;
        out.resize( prev.size() << 1 );

        int offset = 0;
        for ( int prefix = 0; prefix != 2; ++prefix )
        {
            for ( size_t i = 0; i != prev.size(); ++i )
            {
                out[offset+i].resize( prev[i].size() + 1 );
                out[offset+i][0] = prefix;
                std::copy( prev[i].begin(), prev[i].end(), out[offset+i].begin()+1 );
            }
            offset += prev.size();
        }

        return out;
    } // func powerSet

    template <class PrimitivesT> inline typename SimAnnOptProblem<PrimitivesT>::Scalar
    SimAnnOptProblem<PrimitivesT>::getExhaustiveBest( MaskType &opt_mask, std::vector<Scalar> const& desired_angles ) const
    {
        using std::vector;
        float min_e = FLT_MAX;

        vector<MaskType> masks = powerSet( _lines.size() );
        for ( size_t mask_id = 0; mask_id != masks.size(); ++mask_id )
        {
            MaskType &mask = masks[mask_id];
            float e = calcEnergy( mask, desired_angles, NULL, false );
            if ( e < min_e )
            {
                min_e   = e;
                opt_mask = mask;
            }
        }

        return min_e;
    }
} // ... ns am

#endif // __GF2_SIMANNOPTPROBLEM_HPP__
