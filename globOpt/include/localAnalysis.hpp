#ifndef __GF2_LOCALANALYSIS_HPP__
#define __GF2_LOCALANALYSIS_HPP__

#include <fstream>
#include "pcltools/util.hpp"
#include "primitives/houghLine.h"

#include "globFit2.h"

namespace am
{
    template <typename Scalar> inline Eigen::Matrix<Scalar,3,1>
    POS( Eigen::Matrix<Scalar,4,1> const& plane ) { return plane.template segment<3>(0) * (-1.f * plane(3)); }
    template <typename Scalar> inline Eigen::Matrix<Scalar,3,1>
    POS( Eigen::Matrix<Scalar,6,1> const& line ) { return line.template segment<3>(0); }
    template <typename Scalar> inline Eigen::Matrix<Scalar,3,1> DIR( Eigen::Matrix<Scalar,4,1> const& plane ) { return plane.template segment<3>(0); }
    template <typename Scalar> inline Eigen::Matrix<Scalar,3,1>
    DIR( Eigen::Matrix<Scalar,6,1> const& line ) { return line.template segment<3>(3); }

    template <typename Scalar> inline std::pair<Scalar,Scalar> // <score,angle>
    LocalAnalysis::matching( Eigen::Matrix<Scalar,3,1> const& p1, Eigen::Matrix<Scalar,3,1> const& p2, std::vector<Scalar> const& desired_angles )
    {
        Scalar angle     = acos( p1.dot(p2) );
        if ( angle != angle ) angle = static_cast<Scalar>(0);
        while ( angle < 0    ) { std::cerr << "[" << __func__ << "]: " << "angle < 0...adding M_PI\n"; angle += M_PI; }
        while ( angle > M_PI ) { std::cerr << "[" << __func__ << "]: " << "angle > M_PI...subtracting M_PI\n"; angle -= M_PI; }

        Scalar tmp_score = static_cast<Scalar>( 0 );
        std::pair<Scalar,Scalar> min_score( std::numeric_limits<Scalar>::max(), static_cast<Scalar>(0) );
        for ( size_t ang_id = 0; ang_id != desired_angles.size(); ++ang_id )
        {
            // skip parallel
            //if ( desired_angles[ang_id] == static_cast<Scalar>(0) ) continue;

            if ( (tmp_score = fabs(desired_angles[ang_id]-angle)) < min_score.first )
            {
                min_score.first  = tmp_score;
                min_score.second = desired_angles[ang_id];
            }
        }

        return min_score;
    }

    template <typename Scalar> inline Scalar
    angularDifference( Eigen::Matrix<Scalar,6,1> const& l1, Eigen::Matrix<Scalar,6,1> const& l2 )
    {
        // TODO: generalie to desired angles
        Scalar ang = acos( fabs(l1.template segment<3>(3).dot(l2.template segment<3>(3))) );
        Scalar diff = std::min( ang, static_cast<Scalar>(M_PI) - ang);

        // needed!
        if ( diff != diff   ) { diff = static_cast<Scalar>(0); /*std::cerr << "[" << __func__ << "]: " << " correcting nan = acos(" << l1.template segment<3>(3).dot(l2.template segment<3>(3)) << ")\n"; */ }
        // not needed:
        if ( diff <  0      ) { std::cerr << "[" << __func__ << "]: " << " diff < 0\n";       diff = fabs(diff); }
        if ( diff >  M_PI_2 ) { std::cerr << "[" << __func__ << "]: " << " diff > 90\n";      diff = M_PI - diff; }
        if ( diff >  M_PI   ) { std::cerr << "[" << __func__ << "]: " << " diff > 180\n";     diff -= M_PI; }

        return diff;
    }

    template <typename Scalar> inline Scalar
    angularDifference( Eigen::Matrix<Scalar,4,1> const& p1, Eigen::Matrix<Scalar,4,1> const& p2 )
    {
        // TODO: acos
        Scalar ang = acos(fabs(p1.template segment<3>(0).dot(p2.template segment<3>(0))));
        Scalar diff = std::min( ang, static_cast<Scalar>(M_PI) - ang);
        if ( diff != diff ) { /*std::cerr << "[" << __func__ << "]: " << " correcting nan\n";*/ diff = static_cast<Scalar>(0); }
        if ( diff < 0 )     { std::cerr << "[" << __func__ << "]: " << " diff < 0\n"; diff = fabs(diff); }
        if ( diff > M_PI_2 ) { std::cerr << "[" << __func__ << "]: " << " diff > 90\n"; diff = M_PI - diff; }
        if ( diff > M_PI )  { std::cerr << "[" << __func__ << "]: " << " diff > 180\n"; diff -= M_PI; }

        return diff;
    }

    template <class LinesT, class PointsT> inline int
    get_close_points( std::vector<int>        & point_ids
                      , std::vector<int> const& line_ids
                      , LinesT           const& lines
                      , PointsT                 cloud
                      , float            const  scale
                      )
    {
        typedef typename LinesT::value_type::Scalar Scalar;

        // collect point ids for each line in the cluster
        point_ids.reserve( cloud->size() );
        for ( int pid = 0; pid != cloud->size(); ++pid )
        {
            for ( int lid_id = 0; lid_id != line_ids.size(); ++lid_id )
            {
                //Scalar dist = point2line( lines[ line_ids[lid_id] ](), cloud->at(pid).getVector3fMap() );
                Scalar dist = lines[ line_ids[lid_id] ].point3Distance( cloud->at(pid).getVector3fMap() );
                if ( dist < scale )
                {
                    point_ids.push_back( pid );
                    //if ( lines_points ) lines_points
                    break;
                }
            }
        }

        return EXIT_SUCCESS;
    }



    template <  class    LinesT                                        // vector<Primitive>
              , class    PointsT                                       // PointCloud::Ptr
              , typename Scalar = typename LinesT::value_type::Scalar  // float
              , int      Dim = LinesT::value_type::Dim               > // 4 (plane) or 6 (line)
              inline int
    cluster_lines( std::vector<std::vector<int> >        & clusters_lines
                   , std::vector<int>                    & lines_clusters
                   , LinesT                         const& lines
                   , PointsT                               cloud
                   , float                          const  similarity_threshold
                   , Scalar                                (*distance)(Eigen::Matrix<Scalar,Dim,1> const& l1, Eigen::Matrix<Scalar,Dim,1> const& l2)
                   )
    {
        //typedef typename TLines::value_type::Scalar Scalar;

        // init clustering
                  lines_clusters.resize( lines.size(), -1 );
                  clusters_lines.clear();

                  if ( similarity_threshold < 0.f )
                  {
                      for ( int lid = 0; lid != lines.size(); ++lid )
                      {
                          clusters_lines.push_back( { lid } );
                          lines_clusters[ lid ] = clusters_lines.size() - 1;
                      }

                      return EXIT_SUCCESS;
                  }


        // add first line to first cluster
        clusters_lines.push_back( {0} );
        lines_clusters[0] = 0;

        // for each line
        Eigen::Matrix<Scalar,Dim,1> line;
        for ( int lid = 1; lid != lines.size(); ++lid )
        {
            // cache
            line = lines[lid]();

            // calc closest cluster
            Scalar closest_cluster_dist = FLT_MAX;
            int    closest_cluster_id   = -1;
            for ( int cid = 0; cid != clusters_lines.size(); ++cid )
            {
                // calculate line-cluster difference
                Scalar max_diff = -1.f;
                for ( int l = 0; l != clusters_lines[cid].size(); ++l )
                    max_diff = std::max( (*distance)(lines[ clusters_lines[cid][l] ](), line), max_diff );

                // save if cluster closer
                if ( max_diff < closest_cluster_dist )
                {
                    closest_cluster_dist = max_diff;
                    closest_cluster_id   = cid;
                }
            }
            if ( closest_cluster_id < 0 ) { std::cerr << "[" << __func__ << "]: " << "no cluster chosen...aaaaa..." << std::endl; return EXIT_FAILURE; }

            if ( fabs(closest_cluster_dist) > similarity_threshold )
            {
                // new cluster
                clusters_lines.push_back( { lid } );
                lines_clusters[ lid ] = clusters_lines.size()-1;
            }
            else
            {
                // add to similar cluster
                clusters_lines[closest_cluster_id].push_back( lid );
                lines_clusters[ lid ] = closest_cluster_id;
            }
        }

        return EXIT_SUCCESS;
    }

    template <class LinesT, typename Scalar = typename LinesT::value_type::Scalar> inline int
    calculateRepresentatives( LinesT                                & representative_lines
                              , std::vector< std::vector<int> >     & clusters_lines_out
                              , std::vector<int>                    & lines_clusters_out
                              , std::vector<int>                    & cluster_representatives
                              , std::vector<std::vector<int> > const& clusters_lines
                              , LinesT                         const& lines )
    {
        typedef GF2::EigenLine              EigenLine;
        typedef typename LinesT::value_type PrimitiveT;


        cluster_representatives.resize( clusters_lines.size() );
        for ( size_t cid = 0; cid != clusters_lines.size(); ++cid )
        {
            // aggregate
            Eigen::Matrix<Scalar,3,1> representative_dir;
            representative_dir.setZero();
            for ( size_t lid_id = 0; lid_id != clusters_lines[cid].size(); ++lid_id )
            {
                const int lid = clusters_lines[cid][lid_id];
                if ( lines[lid].dir().dot( lines[clusters_lines[cid][0]].dir() ) < 0.f )
                {
                    std::cerr << "lines["<<lid<<"lid].dir(): " << lines[lid].dir().transpose()
                              << " vs. " << lines[clusters_lines[cid][0]].dir().transpose()
                            << ", dot: " << lines[lid].dir().dot( lines[clusters_lines[cid][0]].dir() ) << std::endl;
                    representative_dir += -1.f * lines[lid].dir();
                }
                else
                    representative_dir += lines[lid].dir();
            }
            representative_dir.normalize();

            // replicate angle at each member
            representative_lines.reserve( representative_lines.size() + clusters_lines[cid].size() );
            clusters_lines_out.resize( cid+1 );
            for ( size_t lid_id = 0; lid_id != clusters_lines[cid].size(); ++lid_id )
            {
                const int lid = clusters_lines[cid][lid_id];
                representative_lines.emplace_back( PrimitiveT(POS(lines[lid]()),representative_dir) );

                // save representative line
                if ( lid_id == 0 ) cluster_representatives[cid] = representative_lines.size() - 1;

                // enroll
                clusters_lines_out[cid].push_back( representative_lines.size() - 1 );
                lines_clusters_out     .push_back( cid                             );
            }
        }

        return EXIT_SUCCESS;
    }

    template <class PrimitivesT, typename Scalar>
    inline int
    LocalAnalysis::filter( PrimitiveClustering<PrimitivesT>        & out
                           , PrimitiveClustering<PrimitivesT> const& in
                           , Scalar                           const  working_scale
                           , Scalar                           const  ang_similarity_threshold )
    {
        typedef typename PrimitivesT::value_type PrimitiveT;

        // copy cluster reps first
        out.cluster_representatives.resize( in.cluster_representatives.size(), -1 );
        for ( int cid = 0; cid != in.cluster_representatives.size(); ++cid )
        {
            out.add( in.getClusterRepresentative(cid), cid );
            out.cluster_representatives[cid] = out.lines.size() - 1;
            for ( int d = 0; d != out.lines.back()().rows(); ++d )
                if ( out.lines.back()()(d) != out.lines.back()()(d) )
                    std::cerr << "[" << __func__ << "]: " << "NAN: " << out.lines.back()().transpose() << std::endl;
        }

        // add primitive one after the other
        for ( int lid = 0; lid != in.lines.size(); ++lid )
        {
            // primitive to add
            PrimitiveT const& prim = in.lines[ lid ];

            bool _isnan = false;
            for ( int d = 0; d != prim().rows() && !_isnan; ++d )
                if ( prim()(d) != prim()(d) )
                    _isnan = true;

            if ( _isnan ) continue;

            // check all previously added ones
            bool different = true; size_t lid1 = 0;
            for ( ; lid1 != out.lines.size() && different; ++lid1 )
            {
                PrimitiveT const& other = out.lines[ lid1 ];

                different &= PrimitiveT::different( /*  primitive: */ prim
                                                    , /*    other: */ other
                                                    , /* pos_diff: */ working_scale
                                                    , /* ang_diff: */ ang_similarity_threshold );
            }

            // if still different, keep
            if ( different ) out.add( prim, in.getClusterID(lid) );
//            else
//            {
//                std::cout << "lid1: " << lid1-1 << "/ " << out.lines.size() << std::endl; fflush(stdout);
//                std::cout << prim().transpose() << " was not different from " << out.lines[lid1-1]().transpose()
//                          << ", at scale " << working_scale << ", and angle thresh " << ang_similarity_threshold * 180.f / M_PI
//                          << ", posdiff: " << fabs(prim()(3) - out.lines[lid1-1]()(3)) << "<? " << working_scale
//                          << ", angdiff: " << acos( prim.dir().dot(out.lines[lid1-1].dir()) ) << "<? " << ang_similarity_threshold
//                          << std::endl;
//            }
        } // ... for all lines

        // copy other flags
        out.directions_in_cluster = in.directions_in_cluster;
        out.first_c = 0; // TODO: fix this...

        // report
        Scalar filtered_ratio = 1.f - static_cast<Scalar>(out.lines.size()) / static_cast<Scalar>(in.lines.size());
        if ( filtered_ratio > .3f ) std::cerr << "[" << __func__ << "]: " << "filtered out " << filtered_ratio * 100.f << "%" << std::endl;
        else                        std::cout << "[" << __func__ << "]: " << "filtered out " << filtered_ratio * 100.f << "%" << std::endl;

        return EXIT_SUCCESS;
    }

    template <class LinesT, class PointsT, class Scalar> inline int
    LocalAnalysis::runSimpler(
            std::vector<PrimitiveClustering<LinesT> >  & clusterings
            , PointsT                                    cloud
            , Scalar                              const  working_scale
            , Scalar                              const  ang_similarity_threshold
            , int                                 const  max_neighbourhood_size
            , std::vector<Scalar>                 const &desired_angles
            , std::pair<Scalar,Scalar> (*matching)(Eigen::Matrix<Scalar,3,1> const& p1, Eigen::Matrix<Scalar,3,1> const& p2, std::vector<Scalar> const& desired_angles)
            , int                                 const  addNFriends
            , Scalar                              const  filter_coeff
            , LinesT                                   * input_primitives
            , Scalar                              const  trunc_at
            )
    {
        using std::vector;
        using std::pair;
        using std::deque;
        using std::set;
        using std::cout;
        using std::endl;
        typedef typename LinesT::value_type LineT;
        const int DIM = ((LineT::Dim == 6) ? 2 : 3); // 2D or 3D ?
        //typedef typename LineT::Scalar                      Scalar;
        //typedef typename PointsT::element_type::PointType   PointT;

        // init output
        clusterings.push_back( PrimitiveClustering<LinesT>() );
        clusterings.push_back( PrimitiveClustering<LinesT>() ); // prepare for clusterings[NCID] coming later...
        const int GEN0 = clusterings.size()-2;
        const int GEN1 = clusterings.size()-1;

        // PROPOSE
        if ( !input_primitives )
            LocalAnalysis::propose( clusterings[GEN0].lines, cloud, NULL, max_neighbourhood_size, working_scale );
        else
            clusterings[GEN0].lines = *input_primitives;

        // CLUSTERs[ GEN0 ]
        {
            cluster_lines( clusterings[GEN0].clusters_lines
                           , clusterings[GEN0].lines_clusters
                           , clusterings[GEN0].lines
                           , cloud
                           , ang_similarity_threshold
                           , &angularDifference
                           );

            // they are not even parallel
            clusterings[GEN0].directions_in_cluster = -1;

            std::cout << "clustercount : " << clusterings[GEN0].clusters_lines.size()
                      << ", lines count: " << clusterings[GEN0].lines.size()
                      << ", point_count: " << cloud->size()
                      << "\n";
        }

        // find cluster centre-s
        calculateRepresentatives( clusterings[GEN1].lines
                                  , clusterings[GEN1].clusters_lines
                                  , clusterings[GEN1].lines_clusters
                                  , clusterings[GEN1].cluster_representatives
                                  , clusterings[GEN0].clusters_lines
                                  , clusterings[GEN0].lines );
        clusterings[GEN1].directions_in_cluster = 1;

        // assign points
        vector< vector<int> > clusters_point_ids( clusterings[GEN1].clusters_lines.size() );
        {
            for ( int cid = 0; cid != clusterings[GEN1].clusters_lines.size(); ++cid )
            {
                // call points to clusters
                get_close_points( clusters_point_ids[ cid ]
                                  , /* line_ids: */ clusterings[GEN1].clusters_lines[cid]
                                  , /*    lines: */ clusterings[GEN1].lines
                                  , /*   points: */ cloud
                                  , /*    scale: */ working_scale
                                  );
                std::cout << "cluster " << cid << " has " << clusters_point_ids[cid].size() << " points.\n";
            }
        }


        // replicate lines at close points
#       warning "don't forget this false"
        if ( false && !input_primitives )
        {
            for ( int cid = 0; cid != clusterings[GEN1].clusters_lines.size(); ++cid )
            {
                //float prob = 1000.f / static_cast<float>(clusterings[GEN1].clusters_lines[cid].size()) / clusters_point_ids[cid].size();
                Eigen::Matrix<Scalar,3,1> dir = clusterings[GEN1].getClusterRepresentative( cid ).dir();
                if ( dir(0) != dir(0) || dir(1) != dir(1) || dir(2) != dir(2) )
                        std::cerr << "nan dir..." << std::endl;
                for ( size_t pid_id = 0; pid_id != clusters_point_ids[cid].size(); ++pid_id )
                {
#                   warning "don't forget  this filter"
                    if ( (rand()/static_cast<float>(RAND_MAX)) < 0.9f ) continue;

                    const int pid = clusters_point_ids[cid][pid_id];

                    clusterings[GEN1].lines.emplace_back( LineT(cloud->at(pid).getVector3fMap(),dir) );

                    // enroll
                    clusterings[GEN1].clusters_lines[cid].push_back( clusterings[GEN1].lines.size() - 1 );
                    clusterings[GEN1].lines_clusters     .push_back( cid                                );
                }
            }

            // add local lines to clusters; no new clusters created, only lines added at the back
            {
                // save position of consensus directed lines
                clusterings[GEN0].first_c = clusterings[GEN0].lines.size();
                // reserve
                clusterings[GEN0].lines         .reserve( clusterings[GEN0].lines         .size() + clusterings[GEN1].lines.size() );
                clusterings[GEN0].lines_clusters.reserve( clusterings[GEN0].lines_clusters.size() + clusterings[GEN1].lines.size() );

                if ( clusterings[GEN0].cluster_representatives.size() > 0 ) std::cerr << "[" << __func__ << "]: " << "we don't want this" << std::endl;
                // reserve space for the representatives
                clusterings[GEN0].cluster_representatives.resize( clusterings[GEN1].clusters_lines.size(), -1 );

                // add all consensus lines from clusterings[GEN1]
                for ( size_t lid = 0; lid != clusterings[GEN1].lines.size(); ++lid )
                {
                    // cluster_id
                    const int cid = clusterings[GEN1].lines_clusters[lid];

                    clusterings[GEN0].lines                .push_back( clusterings[GEN1].lines[lid]       ); // copy line
                    clusterings[GEN0].lines_clusters       .push_back( cid                                ); // copy line_id    -> cluster_id
                    clusterings[GEN0].clusters_lines[ cid ].push_back( clusterings[GEN0].lines.size() - 1 ); // copy cluster_id -> line_id
                    // copy representatives, if empty
                    if ( clusterings[GEN0].cluster_representatives[cid] < 0 ) clusterings[GEN0].cluster_representatives[cid] = clusterings[GEN0].lines.size()-1;
                }

                // unsure is this is useful
                clusterings[GEN0].directions_in_cluster = 1;
            }
        }

        // add +-90 to clusters
        if ( addNFriends != 0 )
        {
            std::vector<Scalar> stripped_desired_angles; // remove 0, 180, 360, etc.
            Scalar min_des_ang_spacing = M_PI;
            {
                Scalar tmp;
                for ( size_t i = 0; i != desired_angles.size(); ++i )
                    //if ( desired_angles[i] / M_PI > 1e-3f )
                    {
                        for ( size_t ang_id = 0; ang_id != stripped_desired_angles.size(); ++ang_id )
                            if ( (tmp=fabs(stripped_desired_angles[ang_id] - desired_angles[i])) < min_des_ang_spacing )
                                min_des_ang_spacing = tmp;

                        stripped_desired_angles.push_back( desired_angles[i] );
                    }
            }
            std::cerr << "min_des_ang_spacing: " << min_des_ang_spacing << std::endl;
            std::cout<<"stripped_desired_angles:";for(size_t vi=0;vi!=stripped_desired_angles.size();++vi)std::cout<<stripped_desired_angles[vi]<<" ";std::cout << "\n";


            const int C = clusterings[GEN1].clusters_lines.size();
            vector< deque< std::pair< int, std::pair<Scalar,Scalar>> > > clusters_friends( C ); // [cluster_id] -> { pair<friend_id,<friend_score,best_angle>>, ... }
            {
                const Scalar ang_limit = std::min(min_des_ang_spacing/2.f, trunc_at );
                std::cerr << "[" << __func__ << "]: " << "ang_limit: " << ang_limit << std::endl;

                // foreach cluster
                for ( size_t cid = 0; cid != C-1; ++cid )
                {
                    // get cluster direction
                    Eigen::Matrix<Scalar,3,1> c_dir = clusterings[GEN1].lines[ clusterings[GEN1].clusters_lines[cid][0] ].dir();

                    // foreach other cluster
                    for ( size_t cid2 = cid+1; cid2 != C; ++cid2 )
                    {
                        // get prospective friend's direction
                        Eigen::Matrix<Scalar,3,1> c2_dir = clusterings[GEN1].lines[ clusterings[GEN1].clusters_lines[cid2][0] ].dir();

                        // score relationship (minimum score ~ score to closest angle is returned)
                        pair<Scalar,Scalar> score = matching( c_dir, c2_dir, stripped_desired_angles );
                        if ( score.first > ang_limit )
                        {
                            //continue;
                            if ( input_primitives )
                            {
                                std::cout << "score.first: " << score.first << " > " << ang_limit << " so would skip cid" << cid << "&cid" << cid2 << std::endl;
                            }
                            continue;
                        }

                        /// insert cid2 into cid's friends
                        {
                            // find place
                            typename deque<pair<int,pair<Scalar,Scalar> > >::iterator it = clusters_friends[cid].begin(); // start
                            while ( (it != clusters_friends[cid].end()) && (it->second.first < score.first) )    ++it;  // step until finding something and too big
                            // insert at position
                            clusters_friends[cid].insert( it, pair<int,pair<Scalar,Scalar>>(cid2,score) );
                            // ensure size
                            if ( addNFriends > 0 )
                                while ( clusters_friends[cid].size() > addNFriends ) clusters_friends[cid].pop_back();
                        }

                        /// insert cid into cid2's friends
                        {
                            // find place
                            typename deque<pair<int,pair<Scalar,Scalar> > >::iterator it = clusters_friends[cid2].begin(); // start
                            while ( (it != clusters_friends[cid2].end()) && (it->second.first < score.first) )    ++it;  // step until finding something and too big
                            // insert at position
                            clusters_friends[cid2].insert( it, pair<int,pair<Scalar,Scalar>>(cid,score) );
                            // ensure size
                            if ( addNFriends > 0 )
                                while ( clusters_friends[cid2].size() > addNFriends ) clusters_friends[cid2].pop_back();
                        }
                    }

                    // add self so that c_j direction replicates at centroid later
                    clusters_friends[cid].insert( clusters_friends[cid].end(), pair<int,pair<Scalar,Scalar>>(cid,pair<Scalar,Scalar>(0.f,0.f)) );
                }
                // cid doesn't go to the last one, so add
                clusters_friends.back().insert( clusters_friends.back().end(), pair<int,pair<Scalar,Scalar>>(clusters_friends.size()-1,pair<Scalar,Scalar>(0.f,0.f)) );
            }
#if 0
            // assign points for second run
            {
                // TODO: assign all points to at least one cluster
                if ( input_primitives )
                {
                    std::vector<int> points_lines;
                    am::GlobFit2::assignPoints( points_lines, NULL, clusterings[GEN0].lines, cloud, static_cast<std::vector<Scalar>*>(NULL) );
                    std::cout << "resizing clusters_point_ids from " << clusters_point_ids.size();
                    clusters_point_ids.clear(); clusters_point_ids.resize( clusters_friends.size() );
                    std::cout << " to " << clusters_point_ids.size() << std::endl;

                    int cid, lid;
                    for ( int pid = 0; pid != cloud->size(); ++pid )
                    {
                        lid = points_lines[ pid ];
                        cid = clusterings[GEN0].lines_clusters[ lid ];
                        clusters_point_ids[cid].push_back( pid );
                    }
                }
            }
#endif

            // replicate from GEN1 to GEN0
            vector< vector<int> > gen1_primitives_points;
            am::GlobFit2::assignPoints( /*        points_lines: */ static_cast< std::vector<int>* >(NULL)
                                        , /*              mask: */ static_cast< std::vector<int>* >(NULL)
                                        , /*        primitives: */ clusterings[GEN1].lines
                                        , /*            points: */ cloud
                                        , /* primitives_points: */ &gen1_primitives_points
                                        , /*       dist_thresh: */ working_scale
                                        , /*         sqr_dists: */ static_cast<std::vector<Scalar>* >(NULL) );

            // for each cluster
            for ( size_t cid = 0; cid != clusters_friends.size(); ++cid )
            {
                // for each line - replicate friend clusters direction on the line centroid
                for ( size_t lid_id = 0; lid_id != clusterings[GEN1].clusters_lines[cid].size(); ++lid_id )
                {
                    const int lid = clusterings[GEN1].clusters_lines[cid][lid_id];

                    // points centroid
                    Eigen::Matrix<Scalar,3,1> centroid( Eigen::Matrix<Scalar,3,1>::Zero() );
                    {
                        for ( size_t pid_id = 0; pid_id != gen1_primitives_points[ lid ].size(); ++pid_id )
                        {
                            centroid += cloud->at(gen1_primitives_points[lid][pid_id]).getVector3fMap(); // TODO: clusters_point_ids can have outliers from parallel but close lines
                        }
                        if ( gen1_primitives_points[lid].size() )
                            centroid /= gen1_primitives_points[lid].size();
                        else
                        {
                            std::cerr << "[" << __func__ << "]: " << "wowowo...no points for line?" << std::endl;
                            continue;
                        }
                    } // ... centroid

                    // cache direction
                    Eigen::Matrix<Scalar,3,1> c_dir = DIR( clusterings[GEN1].lines[ clusterings[GEN1].clusters_lines[cid][0] ]() ); // GEN1 has only c_j directions

                    // for each friend
                    for ( typename deque<pair<int,pair<Scalar,Scalar> > >::iterator it = clusters_friends[cid].begin(); it != clusters_friends[cid].end(); ++it )
                    {
                        // debug
                        if ( input_primitives && clusters_friends[cid].size() != input_primitives->size() )
                        {
                            std::cerr << "wtf...cluster_friends[" << cid << "].size(): " << clusters_friends[cid].size() << "\n";
                            std::cout<<"cluster_friends["<<cid<<"]:";
                            for(size_t vi=0;vi!=clusters_friends[cid].size();++vi)std::cout<<clusters_friends[cid][vi].first<<" ";std::cout << "\n";
                        }

                        //  cache
                        int    friend_cid = it->first;
                        Scalar angle      = it->second.second;

                        // calculate rotated direction
                        Eigen::Matrix<Scalar,3,1> dir    = DIR( clusterings[GEN1].lines[ clusterings[GEN1].clusters_lines[friend_cid][0] ]() );
                        if ( angle != static_cast<Scalar>(0) )
                        {
                            Eigen::Matrix<Scalar,3,1> up_dir = dir.cross( c_dir ).normalized();
                            dir = Eigen::AngleAxis<Scalar>( angle, up_dir ) * dir; // ROTATE // TODO: rotate by all possible directions ? (+/-?)
                        }

                        // replicate
                        if ( true || input_primitives ) // add at centroid
                        {
                            clusterings[GEN0].add( LineT(centroid, dir), cid, /* origin: */ friend_cid );
                            std::cout << "added line " << clusterings[GEN0].lines.back()().transpose() << std::endl;
                        }
                        else // add at every point
                        {
                            for ( size_t pid_id = 0; pid_id != clusters_point_ids[ cid ].size(); ++pid_id )
                            {
                                int pid = gen1_primitives_points[lid][pid_id];
                                LineT prim = LineT( cloud->at(pid).getVector3fMap(), dir );

                                clusterings[GEN0].lines              .push_back( prim                               );
                                clusterings[GEN0].clusters_lines[cid].push_back( clusterings[GEN0].lines.size() - 1 );
                                clusterings[GEN0].lines_clusters     .push_back( cid                                );
                            } // ...for pid_id
                        } // ... if at centroid or all points
                    } // ... for each friend
                } // ... for each cluster_member_line
            } // ... for each cluster
        } // ... if addNfriends > 0

        // filter
        if ( filter_coeff > static_cast<Scalar>(0) )
        {
            std::cout << "filter_coeff: " << filter_coeff << std::endl;
            PrimitiveClustering<LinesT> tmp_clustering;
            filter( tmp_clustering, clusterings[GEN0]
                    , working_scale            / static_cast<Scalar>(pow(10,DIM-1)) * filter_coeff
                    , ang_similarity_threshold / static_cast<Scalar>(pow(10,DIM-1)) * filter_coeff );
            std::cout << "tmp_clustering.size(): " << tmp_clustering.lines.size() << std::endl;
            clusterings[GEN0] = tmp_clustering;
        }

        // maintain clustering aux variables
        for ( int cid = 0; cid != clusterings[GEN0].clusters_lines.size(); ++cid )
        {
            if ( (clusterings[GEN0].cluster_representatives.size() <= cid) || (clusterings[GEN0].cluster_representatives[cid] < 0) )
            {
                if ( clusterings[GEN0].clusters_lines[cid].size() )
                {
                    clusterings[GEN0].cluster_representatives.resize( cid+1, -1 );
                    clusterings[GEN0].cluster_representatives[cid] = clusterings[GEN0].clusters_lines[cid][0];
                }
                else
                    std::cerr << "clusterings[GEN0].clusters_lines["<<cid<<"].size() == 0, can't set representative" << std::endl;
            }
        }
        clusterings[GEN0].directions_in_cluster = 1;

        for ( size_t lid = 0; lid != clusterings[GEN0].lines.size(); ++lid )
        {
            //std::cout << clusterings[GEN0].lines[lid]().transpose() << std::endl;
            auto tmp = clusterings[GEN0].lines[lid]();
            if ( (tmp.array() != tmp.array()).all() )
                std::cerr << "WTF: " << tmp.transpose() << std::endl;
            for ( int d = 0; d != tmp.rows(); ++d )
                if ( !std::isfinite(tmp(d)) )
                    std::cerr << "WTF: " << tmp.transpose() << std::endl;

        }

        return EXIT_SUCCESS;
    }

#if 0
    template <class TLines, class PointsT>
    inline int
    LocalAnalysis::run( std::vector<PrimitiveClustering<TLines> > &clusterings, PointsT cloud, float scale_arg )
    {
        using std::vector;
        using std::pair;
        using std::set;
        using std::cout; using std::endl;
        typedef typename TLines::value_type                 TLine;
        typedef typename TLine::Scalar                      Scalar;
        typedef typename PointsT::element_type::PointType   PointT;

        // PARAMs
        const float           similarity_threshold = 0.005f;
        const int             K                    = 10;
        const vector<float>   range( 1, scale_arg );
        //for ( int i = 5; i < 6; ++i ) range.push_back( range.back() * 2.f );
        const float scale = range[0];

        // init output
        clusterings.push_back( PrimitiveClustering<TLines>() );
        clusterings.push_back( PrimitiveClustering<TLines>() ); // prepare for clusterings[NCID] coming later...
        const int GEN0 = clusterings.size()-2;
        const int GEN1 = clusterings.size()-1;

        // PROPOSE
        LocalAnalysis::propose( clusterings[GEN0].lines, cloud, NULL, K, scale );

        // SCALEs
        for ( int scale_id = 0; scale_id != range.size(); ++scale_id )
        {
            // CLUSTERs[ GEN0 ]
            {
                cluster_lines( clusterings[GEN0].clusters_lines
                               , clusterings[GEN0].lines_clusters
                               , clusterings[GEN0].lines
                               , cloud
                               , similarity_threshold
                               , &angularDifference
                               );

                // they are not even parallel
                clusterings[GEN0].directions_in_cluster = -1;

                std::cout << "clustercount : " << clusterings[GEN0].clusters_lines.size()
                          << ", lines count: " << clusterings[GEN0].lines.size()
                          << ", point_count: " << cloud->size()
                          << "\n";
            }

            // NEIGHs
            std::vector<std::vector<int  > > neighs;
            std::vector<std::vector<float> > sqr_dists;
            {
                smartgeometry::getNeighbourhoodIndices( neighs
                                                        , cloud
                                                        , NULL
                                                        , &sqr_dists
                                                        , K   // 10
                                                        , scale // 0.01f
                                                        );
            }

            // CLUSTERs[ GEN1 ]
            // for each line cluster
            vector< vector<int> > clusters_point_ids( clusterings[GEN0].clusters_lines.size() );
            for ( int cid = 0; cid != clusterings[GEN0].clusters_lines.size(); ++cid )
            {
                // call points to clusters
                get_close_points( clusters_point_ids[ cid ]
                                  , /* line_ids: */ clusterings[GEN0].clusters_lines[cid]
                                  , /*    lines: */ clusterings[GEN0].lines
                                  , /*   points: */ cloud
                                  , /*    scale: */ scale / 2.f );
                std::cout << "cluster " << cid << " has " << clusters_point_ids[cid].size() << " points.\n";

                // create local neighbourhood for all points in cluster
                typename PointsT::element_type local_cloud;
                for ( int pid_id = 0; pid_id != clusters_point_ids[cid].size(); ++pid_id )
                {
                    int pid = clusters_point_ids[cid][pid_id];

                    if ( neighs[pid].size() < 2 ) continue;

                    PointT pnt2;
                    for ( int nid = 1; nid != neighs[pid].size(); ++nid )
                    {
                        Eigen::Vector3f pnt = cloud->at(neighs[pid][nid]).getVector3fMap() - cloud->at(pid).getVector3fMap();
                        pnt2.x = pnt(0); pnt2.y = pnt(1); pnt2.z = pnt(2);
                        local_cloud.push_back( pnt2 );
                    }
                }

                // fit lines
                {
                    Eigen::Matrix<Scalar,6,1> line2;
                    int err = smartgeometry::geometry::fitLinearPrimitive( /*           output: */ line2
                                                                , /*         points: */ local_cloud
                                                                , /*          scale: */ scale
                                                                , /*        indices: */ NULL
                                                                , /*    refit times: */ 2
                                                                , /* use input line: */ false
                                                                );
                    if ( (err != EXIT_SUCCESS) || !((line2-line2).array() == (line2-line2).array()).all() )
                    {
                        std::cerr << "line nan: " << line2.transpose() << "continue-ing..." << std::endl;
                        continue;
                    }

                    std::cout << "line fit to " << cid << ": " << line2.transpose() << std::endl;
                    // create new cluster
                    clusterings[GEN1].clusters_lines.push_back( {} );

                    // replicate line
                    for ( int pid_id = 0; pid_id < clusters_point_ids[cid].size(); pid_id += 1 )
                    {
                        // cache
                        int pid     = clusters_point_ids[cid][ pid_id ];

                        // copy direction
                        Eigen::Matrix<Scalar,6,1> line_to_add = line2;
                        // move to point
                        line_to_add.template head<3>() = cloud->at( pid ).getVector3fMap();

                        // save
                        {
                            clusterings[GEN1].lines.push_back( line_to_add );
                            int new_id    = clusterings[GEN1].lines.size()-1;
                            int new_label = clusterings[GEN1].clusters_lines.size()-1;
                            clusterings[GEN1].lines_clusters.push_back( new_label );
                            clusterings[GEN1].clusters_lines.back().emplace_back( new_id );
                        }
                    }
                }

                // they are all parallel
                clusterings[GEN1].directions_in_cluster = 1;
            } // ... clusters gen1

            // fit perpendiculars
            {
                // select best perpendicular cluster dir for each cluster dir
                set<pair<int,int> > perpendiculars;
                set<pair<int,int> > perpendiculars_ordered;
                {
                    const int gen1_cluster_count = clusterings[GEN1].clusters_lines.size();
                    for ( int cid0 = 0; cid0 != gen1_cluster_count; ++cid0 )
                    {
                        // any line from cluster has the same direction...use the first
                        int line0_id = clusterings[GEN1].clusters_lines[cid0].front();

                        int   max_id  = 0;
                        float max_ang = 0.f;
                        for ( int cid1 = 0; cid1 != gen1_cluster_count; ++cid1 )
                        {
                            if ( cid0 == cid1 ) continue;

                            // any line from cluster has the same direction...use the first
                            int   line1_id = clusterings[GEN1].clusters_lines[cid1].front();
                            float ang      = fabs( angularDifference(clusterings[GEN1].lines[line0_id](),
                                                                     clusterings[GEN1].lines[line1_id]()) ); // close to 0 means perpendicular
                            if ( ang > max_ang )
                            {
                                max_ang = ang;
                                max_id  = cid1;
                            }
                        } // ... cid1

                        // create ordered pair
                        pair<int,int> ordered;
                        if ( cid0 < max_id ) ordered = pair<int,int>(  cid0,max_id);
                        else                 ordered = pair<int,int>(max_id,cid0  );

                        // insert pair, if not already there
                        if ( perpendiculars_ordered.find(ordered) == perpendiculars_ordered.end() )
                        {
                            // output
                            perpendiculars.insert( pair<int,int>(  cid0,max_id) );
                            // for duplicate checking
                            perpendiculars_ordered.insert( ordered );
                        }
                    } // ... cid0
                } // ... perpendiculars

                // plot
                for ( auto it = perpendiculars.begin(); it != perpendiculars.end(); ++it  )
                    std::cout << it->first << "-" << it->second << std::endl;

                // CLUSTERs[ GEN2 ] -> create new generation
                clusterings.push_back( clusterings[GEN1] );
                const int GEN2 = clusterings.size()-1;

                // gather local cloud by rotating one of the clusters
                for ( auto it = perpendiculars.begin(); it != perpendiculars.end(); ++it  )
                {
                    // create local neighbourhood for all points in cluster
                    typename PointsT::element_type local_cloud;

                    // get rotation angle in radians
                    //Scalar anglediff = acos( clusterings[GEN1].lines[ clusterings[GEN1].clusters_lines[it->first ][0] ]().template segment<3>(3).dot(
                    //                         clusterings[GEN1].lines[ clusterings[GEN1].clusters_lines[it->second][0] ]().template segment<3>(3) ));
                    auto rot = Eigen::AngleAxis<Scalar>( M_PI_2, Eigen::Vector3f::UnitZ() );

                    // for both clusters in pair
                    for ( int i = 0; i != 2; ++i ) // TODO: do this the other way as well
                    {
                        // cluster id
                        const int cid = i ? it->second : it->first;
                        // for all points in cluster
                        for ( int pid_id = 0; pid_id != clusters_point_ids[cid].size(); ++pid_id )
                        {
                            // get actual point id
                            int pid = clusters_point_ids[cid][pid_id];

                            // skip if no neighbours
                            if ( neighs[pid].size() < 2 ) continue;

                            // add all neighbours to local cloud
                            for ( int nid = 1; nid != neighs[pid].size(); ++nid )
                            {
                                Eigen::Vector3f pnt = cloud->at(neighs[pid][nid]).getVector3fMap() - cloud->at(pid).getVector3fMap();
                                if ( i ) pnt = rot * pnt;

                                PointT pnt2; pnt2.x = pnt(0); pnt2.y = pnt(1); pnt2.z = pnt(2);
                                if ( i ) { pnt2.r = 255.f; }
                                local_cloud.push_back( pnt2 );
                            }
                        }
                    }

                    // fit line
                    {
                        Eigen::Matrix<Scalar,6,1> line2;
                        int err = smartgeometry::geometry::fitLinearPrimitive( /*           output: */ line2
                                                                    , /*         points: */ local_cloud
                                                                    , /*          scale: */ scale
                                                                    , /*        indices: */ NULL
                                                                    , /*    refit times: */ 2
                                                                    , /* use input line: */ false
                                                                    );
                        if ( (err != EXIT_SUCCESS) || !((line2-line2).array() == (line2-line2).array()).all() )
                        {
                            std::cerr << "line nan: " << line2.transpose() << "continue-ing..." << std::endl;
                            continue;
                        }

                        std::cout << "line fit to " << it->first << "-" << it->second << ": " << line2.transpose() << std::endl;
                        // create perpendicular
                        Eigen::Matrix<Scalar,6,1> line3 = line2;
                        line3.template segment<3>(3) = Eigen::AngleAxis<Scalar>(M_PI_2,Eigen::Vector3f::UnitZ()) * line3.template segment<3>(3);

#if 0
                        // debug vis
                        {
                            char tit[255]; sprintf( tit, "%d-%d_0", it->first, it->second );
                            pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer(tit) );
                            vptr->setBackgroundColor( .5, .5, .6 );
                            vptr->addPointCloud( local_cloud.makeShared() );
                            vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3.0 );

                            Eigen::Vector3f start_pt = line2.template segment<3>(3) * -.1f;
                            Eigen::Vector3f end_pt   = line2.template segment<3>(3) * .1f;

                            vptr->addLine( pcl::PointXYZ(start_pt(0),start_pt(1),start_pt(2))
                                           , pcl::PointXYZ(end_pt(0),end_pt(1),end_pt(2))
                                           , 1., 0., 1.
                                           , tit );

                            start_pt = line3.template segment<3>(3) * -.1f;
                            end_pt   = line3.template segment<3>(3) * .1f;

                            sprintf( tit, "%d-%d_1", it->first, it->second );

                            vptr->addLine( pcl::PointXYZ(start_pt(0),start_pt(1),start_pt(2))
                                           , pcl::PointXYZ(end_pt(0),end_pt(1),end_pt(2))
                                           , 0., 0., 1.
                                           , tit );

                            vptr->spin();
                        }
#endif

                        // overwrite directions of first cluster
                        vector<int> &line_ids_0 = clusterings[GEN1].clusters_lines[it->first],
                                    &line_ids_1 = clusterings[GEN1].clusters_lines[it->second];

                        // first.dir=line2.dir, add lines Line(second.pos,line3.dir)
                        for ( int lid_id = 0; lid_id != line_ids_0.size(); ++lid_id )
                        {
                            clusterings[GEN2].lines[ line_ids_0[lid_id] ]().template segment<3>(3) = line2.template segment<3>(3);

                            // copy line position
                            Eigen::Matrix<Scalar,6,1> new_line = clusterings[GEN1].lines[ line_ids_0[lid_id] ]();
                            // overwrite line direction to perpendicular to line3
                            new_line.template segment<3>(3) = line2.template segment<3>(3);
                            // add line
                            clusterings[GEN2].lines.emplace_back( TLine(new_line) );
                            // add to cluster
                            clusterings[GEN2].clusters_lines[it->second].push_back( clusterings[GEN2].lines.size()-1 );
                            clusterings[GEN2].lines_clusters.push_back( it->second );
                        }

                        // second.dir=line3.dir, add lines Line(first.pos,line2.dir)
                        for ( int lid_id = 0; lid_id != line_ids_1.size(); ++lid_id )
                        {
                            // copy direction of second cluster
                            clusterings[GEN2].lines[ line_ids_1[lid_id] ]().template segment<3>(3) = line3.template segment<3>(3);

                            // copy line position
                            Eigen::Matrix<Scalar,6,1> new_line = clusterings[GEN1].lines[ line_ids_1[lid_id] ]();
                            // overwrite line direction to perpendicular to line2
                            new_line.template segment<3>(3) = line3.template segment<3>(3);
                            // add line
                            clusterings[GEN2].lines.emplace_back( TLine(new_line) );
                            // add to cluster
                            clusterings[GEN2].clusters_lines[it->first].push_back( clusterings[GEN2].lines.size()-1 );
                            clusterings[GEN2].lines_clusters.push_back( it->first );
                        }
                    } // ... fit line
                } // ... for perpendiculars

                // they are all parallel
                clusterings[GEN2].directions_in_cluster = 2;

                // create third generation, and replicate directions at points
            } // ... fit perpendiculars
        } // ... scale

        return EXIT_SUCCESS;
    }
#endif

    template <typename Scalar>
    Scalar local_pointPrimitiveDistance(Eigen::Matrix<Scalar,3,1> const& pnt, Eigen::Matrix<Scalar,6,1> const& line )
    {
        return (line.template head<3>() - pnt).cross( line.template segment<3>(3) ).norm();
    }

    template <class TLines, class PointsPtrT>
    inline int
    LocalAnalysis::propose( TLines                  & lines
                            , PointsPtrT              cloud
                            , std::vector<int> const* indices
                            , int                     K
                            , float                   radius
                            //, std::vector<int>      * mapping
                            )
    {
        using std::vector;
        typedef typename TLines::value_type         TLine;
        typedef typename TLine::Scalar              Scalar;
        typedef typename PointsPtrT::element_type   PointsT;

        if ( indices ) { std::cerr << __PRETTY_FUNCTION__ << "]: indices must be NULL, not implemented yet..." << std::endl; return EXIT_FAILURE; }

        // get neighbourhoods
        std::vector<std::vector<int   > > neighs;
        std::vector<std::vector<Scalar> > sqr_dists;
        smartgeometry::getNeighbourhoodIndices( neighs
                                                , cloud
                                                , indices
                                                , &sqr_dists
                                                , K // 5
                                                , radius // 0.02f
                                                , /* soft_radius: */ true
                                                );

        // only use, if more then 10 data-points
        if ( std::count_if( neighs.begin(), neighs.end(), [] (vector<int> const& n1) { return n1.size() > 2; } ) < 2 )
        {
            std::cerr << "[" << __func__ << "]: " << "not enough to work with (<2)...change scale " << radius << std::endl;
            return EXIT_SUCCESS;
        }

        // every point proposes primitive[s] using its neighbourhood
        //const float prob = 300.f / neighs.size();
        int skipped = 0;
        for ( size_t pid = 0; pid != neighs.size(); ++pid )
        {
            //if ( (rand()/static_cast<float>(RAND_MAX)) > prob ) continue; // subsample points

            // can't fit a line to 0 or 1 points
            if ( neighs[pid].size() < 2 )
            {
                ++skipped;
                continue;
            }

            Eigen::Matrix<Scalar,TLine::Dim,1> line;
            int err = smartgeometry::geometry::fitLinearPrimitive<PointsT,Scalar,TLine::Dim>( /*           output: */ line
                                                                                              , /*         points: */ *cloud
                                                                                              , /*          scale: */ radius
                                                                                              , /*        indices: */ &(neighs[pid])
                                                                                              , /*    refit times: */ 2
                                                                                              , /* use input line: */ false
                                                                                              );
            // debug
//            {
//                // debug
//                pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
//                {
//                    vptr->setBackgroundColor( .5, .5, .6 );
//                    PointsPtrT tmp_cloud( new PointsT() );
//                    for ( int nid = 0; nid != neighs[pid].size(); ++nid )
//                        tmp_cloud->push_back( cloud->at(neighs[pid][nid]) );
//                    vptr->addPointCloud( tmp_cloud );

//                }
//                char plane_name[255]; sprintf( plane_name, "plane%d", pid );
//                TLine::draw( line, cloud, radius, &(neighs[pid]), vptr, plane_name, .7, .7, .2, 0 );
//                pcl::ModelCoefficients plane_coeffs;
//                plane_coeffs.values.resize(4); std::copy( line.data(), line.data() + 4, plane_coeffs.values.begin() );

//                sprintf( plane_name, "plane2_%d", pid );
//                vptr->addPlane( plane_coeffs, cloud->at(pid).x, cloud->at(pid).y, cloud->at(pid).z, plane_name, 0 );

//                vptr->spin();
//            }

            if ( err == EXIT_SUCCESS )      lines.emplace_back( TLine(line) );
        }
        std::cout << "[" << __func__ << "]: "
                  << skipped << "/" << neighs.size() << ": " << skipped / static_cast<float>(neighs.size()) * 100.f << "% of points did not produce primitives, so the primitive count is:"
                  << lines.size() << " = " << lines.size() / static_cast<float>(neighs.size()) *100.f << "%" << std::endl;

        return EXIT_SUCCESS;
    }

#if GF2_USE_PCL
    template <class TLines, class PointsT>
    inline int
    LocalAnalysis::display( TLines                                  &lines
                            , PointsT                                cloud
                            , typename TLines::value_type::Scalar    scale
                            , PrimitiveClustering<TLines>                    *p_clustering
                            , const int disp_limit
                            )
    {
        typedef typename TLines::value_type PrimitiveT;

        using std::vector;
        using std::pair;

        auto getCloudName    = [](int vpid) { return "big_cloud"; };
        auto getLineName     = [](int line_id, int vpid) { char tit[255]; sprintf( tit, "line%d", line_id ); return std::string(tit); };

        const int dim = ((PrimitiveT::Dim == 6) ? 2 : 4);

        std::vector<Eigen::Vector3f> colours;
        if ( p_clustering )
        {
            std::cout << "p_clustering->clusters_lines.size(): " << p_clustering->clusters_lines.size() << std::endl;
            colours = am::util::nColoursEigen( p_clustering->clusters_lines.size(), 1.f, true );
        }
        else
        {
            std::cerr << "[" << __func__ << "]: " << "no clustering???" << std::endl;
            return EXIT_FAILURE;
        }


        std::vector<pcl::visualization::PCLVisualizer::Ptr> vptrs;

        // sort by size
        struct CidWithSize {
                int cid;
                int size;
                CidWithSize( int cid, int size) : cid(cid), size(size) {}
                bool operator<(CidWithSize const& b) const { return size > b.size; }
        };
        vector<CidWithSize> sorted;
        {
             sorted.reserve( p_clustering->clusters_lines.size() );
            for ( size_t cid = 0; cid != p_clustering->clusters_lines.size(); ++cid )
            {
                sorted.push_back( CidWithSize(cid,p_clustering->clusters_lines[cid].size()) );
            }
            std::sort( sorted.begin(), sorted.end() );
        }

        const int N = p_clustering ? std::min(sorted.size(),static_cast<size_t>(disp_limit)) : 1;
        int dx = 0, dy = 0;
        while ( dx * dy < N )
            if ( dx < dy * 1.6 ) ++dx; else ++dy;
        const int w = 1920/dx, h = 1200/dy;

        for ( int cid_id = 0; cid_id != N; ++cid_id )
        {
            const int cid = sorted[cid_id].cid;

            // ang diff
            float ang_diff = 0.f;
            if ( p_clustering )
            {
                float tmp;
                for ( size_t k = 0; k != p_clustering->clusters_lines[cid].size()-1; ++k )
                    for ( size_t k1 = k+1; k1 != p_clustering->clusters_lines[cid].size(); ++k1 )
                    {
                        tmp = acos( lines[k].dir().dot( lines[k1].dir() ) );
                        if ( tmp != tmp ) tmp = 0.f;

                        if ( tmp > ang_diff )
                            ang_diff = tmp;
                    }

            }

            char tit[255];
            sprintf( tit, "cid%d size: %d, max_ang_diff: %f radians"
                     , cid
                     , static_cast<int>(p_clustering ? p_clustering->clusters_lines[cid].size() : lines.size())
                     , ang_diff );

            pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer(tit) );
            vptrs.push_back(vptr);

            std::string cloud_name = getCloudName( 0 );
            vptr->setBackgroundColor( .6, .6, .6, 0 );
            vptr->setSize( w, h);
            vptr->setPosition( (cid_id % dx) * w, (cid_id / dx) * h );
            vptr->addPointCloud( cloud, cloud_name, 0 );
            vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3., cloud_name, 0 );

            // show lines
            vector<pcl::PointXYZ> endpoints; endpoints.reserve( lines.size()*dim );
            for ( size_t k = 0; k != lines.size(); ++k )
            {
                std::vector<typename PointsT::element_type::PointType> minMax;
                int err = lines[k].getExtent( minMax
                                              , cloud
                                              , scale
                                              //, labels ? &(labels_points[k])
                                              //         : indices ? &(indices->indices) :
                                              , NULL );

                if ( err != EXIT_SUCCESS )  continue;

                double colour[3] = { 0., .3, 1. };
                if ( p_clustering )
                {
                    int cluster_id = p_clustering->lines_clusters[ k ];
                    if ( cluster_id!=cid ) continue;
                    colour[0] = colours[cluster_id](0);
                    colour[1] = colours[cluster_id](1);
                    colour[2] = colours[cluster_id](2);
                }
                std::vector<pcl::PointXYZ> ps(dim);
                if ( dim == 2 )
                {
                    Eigen::Vector3f ep0 = minMax[0].getVector3fMap() - .2f * (minMax[1].getVector3fMap() - minMax[0].getVector3fMap());
                    ps[0] = pcl::PointXYZ( ep0(0), ep0(1), ep0(2) );
                    Eigen::Vector3f ep1 = minMax[1].getVector3fMap() - .2f * (minMax[0].getVector3fMap() - minMax[1].getVector3fMap());
                    ps[1] = pcl::PointXYZ( ep1(0), ep1(1), ep1(2) );
                }
                else if ( dim == 4 )
                {
                    pcl::PointXYZ tmp;
                    for ( int corner_id = 0; corner_id != minMax.size(); ++corner_id )
                    {
                        tmp.x = minMax[corner_id].x; tmp.y = minMax[corner_id].y; tmp.z = minMax[corner_id].z;
                        ps[corner_id] = tmp ;
                    }
                }
                else std::cerr << "[" << __func__ << "]: " << " not ok dim" << std::endl;

                PrimitiveT::draw( ps, vptr, getLineName(k,0), colour[0], colour[1], colour[2], 0 );
//                vptr->addLine( p0
//                               , p1
//                               , colour[0], colour[1], colour[2]
//                               , getLineName( k, 0 )
//                               , 0 );
            }
            vptr->resetCamera();
            vptr->spinOnce();
        }
        vptrs.back()->spin();

        return EXIT_SUCCESS;
    }

#endif

    template <class PointsT>
    inline int
    LocalAnalysis::plot( PointsT cloud )
    {
        using std::vector;

        const int K = 10;
        vector<float> range(1,.001f);
        for ( int i = 0; i < 10; ++i ) range.push_back( range.back() * 2.f );

        std::ofstream gpfile("locals.gp");
        if ( !gpfile.is_open() )
        {
            std::cerr << "[" << __func__ << "]: " << "cannot open locals.gp" << std::endl;
            return EXIT_FAILURE;
        }
        int rids_written = 0;

        std::ofstream f3d( "locals3D.txt" );
        if ( !f3d.is_open() )
        {
            std::cerr << "[" << __func__ << "]: " << "cannot open locals3D.txt" << std::endl;
            return EXIT_FAILURE;
        }

        for ( int rid = 0; rid != range.size(); ++rid )
        {
            float scale = range[rid];
            std::vector<std::vector<int  > > neighs;
            std::vector<std::vector<float> > sqr_dists;
            smartgeometry::getNeighbourhoodIndices( neighs
                                                    , cloud
                                                    , NULL
                                                    , &sqr_dists
                                                    , K   // 10
                                                    , scale // 0.01f
                                                    );
            // only plot, if more then 10 data-points
            if ( std::count_if( neighs.begin(), neighs.end(), [] (vector<int> const& n1) { return n1.size() > 2; } ) < 10 )
                continue;
            //if ( std::max_element( neighs.begin(), neighs.end(), [] (vector<int> const& n1, vector<int> const& n2) { return n1.size() < n2.size(); })->size() < 2 )
            //    continue;

            char fname[255];
            sprintf( fname, "locals_%4.3f.txt", scale );

            std::ofstream f;
            f.open(fname);
            if ( !f.is_open() )
            {
                std::cerr << "[" << __func__ << "]: " << "could not open locals.txt" << std::endl;
                return EXIT_FAILURE;
            }

            int counter = 0;
            for ( int pid = 0; pid != neighs.size(); ++pid )
                //for ( int pid = 0; pid != 10; ++pid )
            {
                //f << plot for [IDX=0:1] 'test.dat' i IDX u 1:2 w lines title columnheader(1)
                if ( neighs[pid].size() < 2 ) continue;

                ++counter;

                f << "\"pid" << pid << "\"" << std::endl;
                for ( int nid = 1; nid != neighs[pid].size(); ++nid )
                {
                    Eigen::Vector3f p = (cloud->at( neighs[pid][nid] ).getVector3fMap() - cloud->at(pid).getVector3fMap());
                    f << pid * K + nid << " "
                      << p.transpose() << "\n";
                    f3d //<< pid * K + nid << " "
                        << p(0) << " " << p(1) << " " << scale
                        << "\n";
                }
                f << "\n\n";
            }

            f.close();
            gpfile << "set term wxt " << rids_written << " persist\n" // position " << 10 + rids_written * 640 << "," << (rids_written / 3) * 490 << "\n"
                   << "unset key;\n"
                   << "plot for [IDX=0:" << counter-1 << "] '" << fname << "' i IDX u 2:3 w points title columnheader(1)\n";

            ++rids_written;

            std::cout << "gnuplotN: " << counter << ", at scale " << scale << std::endl;
        }

        gpfile.close();
        f3d.close();

        return EXIT_SUCCESS;
    }

} // ns am

#endif // __GF2_LOCALANALYSIS_HPP__
