#ifndef __GF2_CANDIDATEGENERATOR_H__
#define __GF2_CANDIDATEGENERATOR_H__

#include <opencv2/core/core.hpp>        // Mat
#include <opencv2/highgui/highgui.hpp>  // imread
#include <opencv2/imgproc/imgproc.hpp>  // cvtColor

#if GF2_USE_PCL
#   include "pcl/common/intersections.h"   // lineWithLineIntersection
#   include "pcl/kdtree/kdtree.h"          // nearestneighboursearch
#   include "globfit2/processing/util.hpp"            // addGaussianNoise, fitLinearPrimitive
#endif // GF2_USE_PCL

#include "globfit2/optimization/patchDistanceFunctors.h" // FullLinkagePointPatchDistanceFunctor
#include "globfit2/parameters.h"        // CandidateGeneratorParams

namespace GF2
{


    // generate
    /// (1) patchify
    //// (1.1) propose
    //// (1.2) agglomerative
    /// (2)
    class CandidateGenerator
    {
        public:
            /*! \brief Main functionality to generate lines from points.
             *
             *  \tparam _PointPrimitiveDistanceFunctorT Concept: \ref MyPointPrimitiveDistanceFunctor.
             *  \tparam _PrimitiveT                     Concept: \ref GF2::LinePrimitive2.
             */
            template <  class       PrimitivePrimitiveAngleFunctorT // concept: energyFunctors.h::PrimitivePrimitiveAngleFunctor
                      , class       _PointPrimitiveDistanceFunctorT
                      , class       _PrimitiveT
                      , class       PrimitiveContainerT             // concept: std::vector<std::vector<LinePrimitive2>>
                      , class       PointContainerT                 // concept: std::vector<PointPrimitive>
                      , typename    Scalar
                      >
            static inline int
            generate( PrimitiveContainerT                     &  out_lines
                    , PrimitiveContainerT              const&  in_lines
                    , PointContainerT                  const&  points // non-const to be able to add group tags
                    , Scalar                           const   scale
                    , std::vector<Scalar>              const&  angles
                    , CandidateGeneratorParams<Scalar> const&  params );
#if GF2_WITH_SAMPLE_INPUT
            //! \brief image_2_2DCloud
            template <  class       PointAllocatorFunctorT
                      , class       PointContainerT>
            static inline int
            image_2_2DCloud( PointContainerT         & cloud
                             , std::string             img_path
                             , int                     N_samples
                             , float const             Z
                             , float const             scale );

            //! \brief image_2_2DCloud
            template <  class PointAllocatorFunctorT
                      , class PointContainerT>
            static inline int
            image_2_2DCloud( PointContainerT         & cloud
                           , cv::Mat         const & img
                           , int   const             N_samples
                           , float const             Z
                           , float const             scale );
#endif // GF2_WITH_SAMPLE_INPUT
    }; //...class CandidateGenerator
} // ...ns::GF2

//_____________________________________________________________________________________________________________________
// HPP
#include "globfit2/my_types.h"          //  PCLPointAllocator
#include "globfit2/processing/util.hpp" // calcPopulations()
#include "globfit2/optimization/segmentation.h" // patchify()

namespace GF2
{
    inline bool equal2D( int l0, int l1, int l2, int l3 ) { return (l0 == l2) && (l1 == l3); }

    /*! \brief Called from \ref CandidateGenerator::generate() via \ref processing::filterPrimitives, decides if a primitive has the GID looked for.
     * \tparam _PrimitiveT Concept: \ref GF2::LinePrimtive2.
     */
    template <class _PointPrimitiveDistanceFunctor, class _PrimitiveT, class _PointPrimitiveT, typename _Scalar>
    struct NearbyPrimitivesFunctor
    {
            NearbyPrimitivesFunctor( _PointPrimitiveT const& point, _Scalar scale ) : _point(point), _scale(scale) {}

            inline int eval( _PrimitiveT const& prim, int const lid )
            {
                // store if inside scale
                if ( _PointPrimitiveDistanceFunctor::template eval<_Scalar>(_point, prim) < _scale )
                {
                    _gidLids.push_back( std::pair<int,int>(prim.getTag(_PrimitiveT::GID), lid) );
                } //...if dist

                return 0;
            }

            _PointPrimitiveT                  _point;
            _Scalar                           _scale;
            std::vector< std::pair<int,int> > _gidLids; //!< Output unique IDs of primitives.
    }; //NearbyPrimitivesFunctor


    /*! \brief Called from getOrphanGids via \ref processing::filterPrimitives, decides if a primitive has the GID looked for.
     */
    template <class _PrimitiveT>
    struct FindGidFunctor
    {
            FindGidFunctor( int gid ) : _gid(gid) {}
            inline int eval( _PrimitiveT const& prim, int /*lid*/ ) const { return prim.getTag(_PrimitiveT::GID) == _gid; }
            int  _gid; //!< \brief GID looked for.
    }; //...FindGidFunctor

    /*! \brief Gathers groupIDs that don't have primitives selected for them.
     *
     *  \tparam _PrimitiveT           Concept: \ref GF2::LinePrimtive2.
     *  \tparam _inner_const_iterator Concept: _PrimitiveContainerT::value_type::const_iterator if vector, _PrimitiveContainerT::mapped_type::const_iterator if map.
     *  \tparam _PidContainerT        Container holding the orphan point ids. Concept: std::set<int>.
     *  \tparam _GidContainerT        Container holding the patch ids (GIDs) with orphan points. Concept: std::set<int>.
     *  \tparam _PointContainerT      Concept: vector<_PointPrimitiveT>.
     *  \tparam _PrimitiveContainerT  Concept: vector< vector< _PrimitiveT> >
     *
     * \param[in] point_gid_tag      Identifies the field, where the GID is stored in the pointPrimitives. Concept: _PointPrimitiveT::GID
     */
    template < class _PrimitiveT
             , class _inner_const_iterator
             , class _PidContainerT
             , class _GidContainerT
             , class _PointContainerT
             , class _PrimitiveContainerT
             >
    inline int getOrphanGids( _GidContainerT            & orphan_gids
                            , _PointContainerT     const& points
                            , _PrimitiveContainerT const& prims
                            , int                  const  point_gid_tag
                            , _PidContainerT            * orphan_pids   = NULL )
    {
        // select unassigned points
        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            int gid = points[pid].getTag( point_gid_tag );

            FindGidFunctor<_PrimitiveT> functor(gid);
            if ( !processing::filterPrimitives< _PrimitiveT, _inner_const_iterator>
                                              ( prims, functor )                                )
            {
                orphan_gids.insert( gid );
                if ( orphan_pids )
                    orphan_pids->insert( pid );
            }
        }
        return EXIT_SUCCESS;
    }


    template <  class       _PrimitivePrimitiveAngleFunctorT
              , class       _PointPrimitiveDistanceFunctorT
              , class       _PrimitiveT
              , class       _PrimitiveContainerT
              , class       _PointContainerT
              , typename    _Scalar> int
    CandidateGenerator::generate( _PrimitiveContainerT                   & out_lines
                                , _PrimitiveContainerT              const& in_lines
                                , _PointContainerT                  const& points   // non-const to be able to add group tags
                                , _Scalar                           const  scale
                                , std::vector<_Scalar>              const& angles
                                , CandidateGeneratorParams<_Scalar> const& params )
    {
        typedef typename _PointContainerT::value_type                     _PointPrimitiveT;
        typedef typename _PrimitiveContainerT::const_iterator             outer_const_iterator;
        //typedef typename outer_const_iterator::value_type::const_iterator inner_const_iterator;
        typedef typename _PrimitiveContainerT::mapped_type::const_iterator inner_const_iterator;

        if ( out_lines.size() ) std::cerr << "[" << __func__ << "]: " << "warning, out_lines not empty!" << std::endl;
        if ( params.patch_population_limit <= 0 ) { std::cerr << "[" << __func__ << "]: " << "error, popfilter is necessary!!!" << std::endl; return EXIT_FAILURE; }

        // (2) Mix and Filter
        // count patch populations
        GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        {
            processing::getPopulations( populations, points );
        }

        // convert local fits
        const _Scalar angle_limit( params.angle_limit / params.angle_limit_div );
        int           nlines = 0;

        // filter already copied directions
        std::map< int, std::set<int> > copied; // [gid] = [ dir0, dir 1, dir2, ... ]

        int gid0, gid1, dir_gid0, dir_gid1;
        gid0 = gid1 = dir_gid0 = dir_gid1 = -1; // group tags cached

        // The following cycles go through each patch and each primitive in each patch and pair them up with each other patch and the primitives within.
        // So, 00 - 00, 00 - 01, 00 - 02, ..., 01 - 01, 01 - 02, 01 - 03, ... , 10 - 10, 10 - 11, 10 - 12, ..., 20 - 20, 20 - 21,  ... (ab - cd), etc.
        // for a = 0 {   for b = 0 {   for c = a {   for d = b {   ...   } } } }

        // OUTER0 (a)
        for ( outer_const_iterator outer_it0  = in_lines.begin();
                                   outer_it0 != in_lines.end();
                                 ++outer_it0 )
        {
            // INNER1 (b)
            gid0 = -2; // -1 is unset, -2 is unread
            for ( inner_const_iterator inner_it0  = containers::valueOf<_PrimitiveT>(outer_it0).begin();
                                       inner_it0 != containers::valueOf<_PrimitiveT>(outer_it0).end();
                                     ++inner_it0 )
            {
                // cache outer primitive
                _PrimitiveT const& prim0 = containers::valueOf<_PrimitiveT>( inner_it0 );
                dir_gid0                 = prim0.getTag( _PrimitiveT::DIR_GID );

                // cache group id of patch at first member
                if ( gid0 == -2 )
                    gid0 = prim0.getTag( _PrimitiveT::GID ); // store gid of first member in patch
                else if ( prim0.getTag(_PrimitiveT::GID) != gid0 )
                    std::cerr << "[" << __func__ << "]: " << "Not good, prims under one gid don't have same GID..." << std::endl;


                // OUTER1 (c)
                gid1 = -2; // -1 is unset, -2 is unread
                for ( outer_const_iterator outer_it1  = in_lines.begin();
                                           outer_it1 != in_lines.end();
                                         ++outer_it1 )
                {
                    if ( outer_it1 == in_lines.begin() )
                         std::advance( outer_it1, std::distance<outer_const_iterator>( in_lines.begin(), outer_it0) );

                    // INNER1 (d)
                    for ( inner_const_iterator inner_it1  = (outer_it0 == outer_it1) ? inner_it0
                                                                                     : containers::valueOf<_PrimitiveT>(outer_it1).begin();
                                               inner_it1 != containers::valueOf<_PrimitiveT>(outer_it1).end();
                                             ++inner_it1 )
                    {
                        // cache inner primitive
                        _PrimitiveT const& prim1 = containers::valueOf<_PrimitiveT>( inner_it1 );
                        dir_gid1                 = prim1.getTag( _PrimitiveT::DIR_GID );

                        // cache group id of patch at first member
                        if ( gid1 == -2 )
                            gid1 = prim1.getTag( _PrimitiveT::GID );

                        // check, if same line (don't need both copies then, since 01-01 =~= 01-01)
                        const bool same_line = ((outer_it0 == outer_it1) && (inner_it0 == inner_it1));

                        bool add0 = same_line, // add0: new primitive at location of prim0, with direction from prim1. We need to keep a copy, so true if same_line.
                             add1 = false;     // add1: new primitive at location of prim1, with direction from prim0

                        // find best rotation id and value
                        int     closest_angle_id = 0;
                        _Scalar closest_angle    = _Scalar( 0 );
                        {
                            _Scalar angdiff    = _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>( prim0, prim1, angles, &closest_angle_id );
                            closest_angle      = angles[ closest_angle_id ];

                            // decide based on rotation difference (RECEIVE_SIMILAR)
                            bool    close_ang  = angdiff < angle_limit;
                            add0              |= close_ang; // prim0 needs to be close to prim1 in angle
                            add1              |= close_ang; // prim1 needs to be close to prim0 in angle
                        }

                        if ( params.small_mode == CandidateGeneratorParams<_Scalar>::SmallPatchesMode::RECEIVE_ALL )
                        {
                            add0 |= (populations[gid0].size() < params.patch_population_limit); // prim0 needs to be small to copy the dir of prim1.
                            add1 |= (populations[gid1].size() < params.patch_population_limit); // prim1 needs to be small to copy the dir of prim0.
                        }
                        else if ( params.small_mode == CandidateGeneratorParams<_Scalar>::SmallPatchesMode::IGNORE )
                        {
                            // both need to be large to work
                            add0 &= add1 &= (populations[gid0].size() >= params.patch_population_limit) && (populations[gid1].size() >= params.patch_population_limit);
                        }

                        // store already copied pairs
                        {
                            if ( copied[gid0].find(dir_gid1) != copied[gid0].end() )
                                add0 = false;

                            if ( copied[gid1].find(dir_gid0) != copied[gid1].end() )
                                add1 = false;
                        }

                        if ( !add0 && !add1 )
                            continue;

                        Eigen::Matrix<_Scalar,3,1> dir0 = prim1.dir(),
                                                   dir1 = prim0.dir();
                        if ( (closest_angle >= _Scalar(0)) && (closest_angle < M_PI) )
                        {
                            dir0 = Eigen::AngleAxisf(-closest_angle, Eigen::Matrix<_Scalar,3,1>::UnitZ() ) * dir0;
                            dir1 = Eigen::AngleAxisf( closest_angle, Eigen::Matrix<_Scalar,3,1>::UnitZ() ) * dir1;
                        }

                        // copy line from pid to pid1
                        if ( add0 )
                        {
                            // prepare
                            _PrimitiveT cand0 = _PrimitiveT( prim0.pos(), dir0 );
                            cand0.setTag( _PrimitiveT::GID    , prim0.getTag(_PrimitiveT::GID) );
                            cand0.setTag( _PrimitiveT::DIR_GID, prim1.getTag(_PrimitiveT::DIR_GID) ); // recently changed this from GID
                            // insert
                            int check = containers::add( out_lines, gid0, cand0 ).getTag( _PrimitiveT::GID );
                            ++nlines;

                            // debug
                            if ( check != cand0.getTag(_PrimitiveT::GID) ) std::cerr << "gid tag copy check failed" << std::endl;
                            int tmp_size = copied[ cand0.getTag(_PrimitiveT::GID) ].size();

                            // keep track of instances
                            copied[ cand0.getTag(_PrimitiveT::GID) ].insert( cand0.getTag(_PrimitiveT::DIR_GID) );

                            // debug
                            if ( copied[cand0.getTag(_PrimitiveT::GID)].size() == tmp_size )
                                std::cerr << "[" << __func__ << "][" << __LINE__ << "]: NOOOOO insertion, should not happen"
                                          << cand0.getTag( _PrimitiveT::GID ) << ", " << cand0.getTag( _PrimitiveT::DIR_GID )
                                          << std::endl;
                        }

                        if ( !same_line && add1 ) // add other copy only, if it's not the same line (we don't want duplicates)
                        {
                            _PrimitiveT cand1 = _PrimitiveT( prim1.pos(), dir1 );
                            cand1.setTag( _PrimitiveT::GID    , prim1.getTag(_PrimitiveT::GID    ) );
                            cand1.setTag( _PrimitiveT::DIR_GID, prim0.getTag(_PrimitiveT::DIR_GID) );

                            // copy line from pid1 to pid
                            int check = containers::add( out_lines, gid1, cand1 ).getTag( _PrimitiveT::DIR_GID );
                            ++nlines;

                            // debug
                            if ( check != cand1.getTag(_PrimitiveT::DIR_GID) ) std::cerr << "dirgid tag copy check failed" << std::endl;
                            int tmp_size = copied[ cand1.getTag(_PrimitiveT::GID) ].size();

                            // keep track of instances
                            copied[ cand1.getTag(_PrimitiveT::GID) ].insert(  cand1.getTag(_PrimitiveT::DIR_GID) );

                            // debug
                            if ( copied[cand1.getTag(_PrimitiveT::GID)].size() == tmp_size )
                                std::cerr << "[" << __func__ << "][" << __LINE__ << "]: NOOOOO insertion, should not happen for " << cand1.getTag( _PrimitiveT::GID ) << ", " << cand1.getTag( _PrimitiveT::DIR_GID ) << std::endl;
                        } //...if add1
                    } //...for l3
                } //...for l2
            } //...for l1
        } //...for l0

        // log
        std::cout << "[" << __func__ << "]: " << "finished generating, we now have " << nlines << " candidates" << std::endl;

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::generate()

#if 0
    int addPrimitivesToSmallPatches()
    {
        // preprocess: add dummy primitives to small clusters so that they can receive others
        if ( params.small_mode == CandidateGeneratorParams<_Scalar>::SmallPatchesMode::RECEIVE_ALL )
        {
            std::set<int> orphan_gids, orphan_pids;
            getOrphanGids<_PrimitiveT, inner_const_iterator>( orphan_gids, points, in_lines, _PointPrimitiveT::GID, &orphan_pids );

            // debug
            std::cout << "orphan_gids: ";
            for ( auto it = orphan_gids.begin(); it != orphan_gids.end(); ++it )
            {
                std::cout << *it << ", ";
            }
            std::cout << std::endl;

            // for all points
            std::map< int, std::set<int> > copied;
            std::set<int>::const_iterator pid_it_end = orphan_pids.end();
            for ( std::set<int>::const_iterator pid_it = orphan_pids.begin(); pid_it != pid_it_end; ++pid_it )
            {
                const int gid = points[*pid_it].getTag( _PointPrimitiveT::GID );

                // get all adopters
                NearbyPrimitivesFunctor<_PointPrimitiveDistanceFunctorT, _PrimitiveT, _PointPrimitiveT, _Scalar> functor( points[*pid_it], params.scale );
                processing::filterPrimitives<_PrimitiveT, inner_const_iterator>( in_lines, functor );

                // debug
                std::cout << "[" << __func__ << "]: " << "point" << *pid_it << "(" << gid << ") has the following adopters: ";
                for ( std::vector<std::pair<int,int> >::const_iterator it = functor._gidLids.begin(); it != functor._gidLids.end(); ++it )
                {
                    _PrimitiveT const& prim( in_lines.at(it->first).at(it->second) );
                    const int dir_gid = prim.getTag(_PrimitiveT::DIR_GID);
                    if ( copied[gid].find( dir_gid ) == copied[gid].end() )
                    {
                        // create candidate
                        _PrimitiveT cand( /* pos: */ processing::getCentroid<_Scalar>(points,populations[gid]) // todo: cache centroids
                                          , /* dir: */ prim.template dir() );
                        // add group ids
                        cand.setTag( _PrimitiveT::GID    , gid     ); // from point patch
                        cand.setTag( _PrimitiveT::DIR_GID, dir_gid ); // from explaining primitive
                        // add to candidates
                        containers::add( out_lines, gid, cand );
                        // remember patch-direction combination
                        copied[gid].insert( dir_gid );
                        ++nlines;

                        // debug
                        std::cout << "!" << prim.getTag(_PrimitiveT::GID) << "," << prim.getTag(_PrimitiveT::DIR_GID) << "!; ";
                    }
                    else
                        std::cout << "(" << prim.getTag(_PrimitiveT::GID) << "," << prim.getTag(_PrimitiveT::DIR_GID) << "); ";
                }
                std::cout << std::endl;


            }
        }
    }
#endif

#if GF2_WITH_SAMPLE_INPUT
    /**
     * @brief GlobFit2::image_2_2DCloud Randomly samples a 2D gray image with 3D points at black pixels to create a 3D cloud.
     * @param cloud
     * @param img_path
     * @param N_samples
     * @param Z
     * @return EXIT_SUCCESS/EXIT_FAILURE
    */
    template <  class PointAllocatorFunctorT
              , class PointContainerT> int
    CandidateGenerator::image_2_2DCloud( PointContainerT              &cloud
                                         , std::string                 img_path
                                         , int                         N_samples
                                         , const float                 Z
                                         , const float scale )
    {
        return image_2_2DCloud( cloud, cv::imread(img_path,cv::IMREAD_UNCHANGED), N_samples, Z, scale );
    } // ...CandidateGenerator::image_2_2DCloud()

    template <  class PointAllocatorFunctorT
              , class PointContainerT> int
    CandidateGenerator::image_2_2DCloud( PointContainerT                &cloud
                                         , cv::Mat                const & img
                                         , const int                      N_samples
                                         , const float                    Z
                                         , const float                    scale )
    {
        typedef typename PointContainerT::PointType PointT;
        //   const float scale = 1.f;

//        if ( !cloud )
//            cloud = pcl::PointCloud<MyPoint>::Ptr( new pcl::PointCloud<MyPoint>() );
        cloud.reserve( N_samples );

        cv::Mat gray;
        if ( img.channels() > 1 )   cv::cvtColor( img, gray, cv::COLOR_RGB2GRAY );
        else                        gray = img;

        cv::threshold( gray, gray, 200., 255., cv::THRESH_BINARY );
        cv::Mat sampled( gray.clone() ); sampled.setTo(0);

        int   pix_count = gray.cols * gray.rows - cv::countNonZero( gray );
        float prob      = N_samples  / (float)pix_count;
        while ( cloud.size() < static_cast<size_t>(N_samples) )
        {
            pix_count = gray.cols * gray.rows - cv::countNonZero( gray );
            prob      = (N_samples - cloud.size())  / (float)pix_count;

            for ( int y = 0; (y < gray.rows) && (cloud.size() != (size_t)N_samples); ++y  )
            {
                for ( int x = 0; (x < gray.cols) && (cloud.size() != (size_t)N_samples); ++x )
                {
//                    float depth = (Eigen::Vector3f() << x, (gray.rows-y), 0).finished().norm();
//                    if ( RANDF() > depth / scale / sqrt(2.f) ) continue;

                    if ( (gray.at<uchar>(y,x) < 255) && (rand()/static_cast<double>(RAND_MAX) < prob) && (sampled.at<uchar>(y,x) == 0) )
                    {
                        //std::cout << "cv::sum( gray(cv::Range(y-1,y+1), cv::Range(x-1, x+1))): " << cv::sum( sampled(cv::Range(std::max(0,y-1),std::min(gray.rows,y+1)),
                        //                                                                                             cv::Range(std::max(0,x-1),std::min(gray.cols,x+1))) )[0] << std::endl;
                        if ( cv::sum( sampled(cv::Range(std::max(0,y-1),std::min(gray.rows,y+1)),
                                              cv::Range(std::max(0,x-1),std::min(gray.cols,x+1))) )[0] >= 3 ) continue;

//                        MyPoint pnt;
//                        pnt.x = x             / (float)gray.cols * scale;
//                        pnt.y = (gray.rows-y) / (float)gray.rows * scale;
//                        pnt.z = Z;

//                        am::util::pcl::setPointColor( pnt, 1, 0, 0 );
//                        cloud.push_back( pnt );


                        cloud.push_back( PointAllocatorFunctorT::template create<PointT>(
                                             (Eigen::Matrix<float,3,1>() <<
                                              x             / (float)gray.cols * scale,
                                              (gray.rows-y) / (float)gray.rows * scale,
                                              Z
                                              ).finished()
                                             )
                                         );

                        sampled.at<uchar>(y,x) = 1;
                    }
                }
            }
        }

        cv::imshow( "img", gray );
        cv::waitKey( 20 );

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::image_2_2DCloud()
#endif // GF2_WITH_SAMPLE_INPUT
} // ... ns GF2


#endif // __GF2_CANDIDATEGENERATOR_H__

