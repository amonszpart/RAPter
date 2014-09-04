#ifndef __GF2_CANDIDATEGENERATOR_H__
#define __GF2_CANDIDATEGENERATOR_H__

#include <opencv2/core/core.hpp>        // Mat
#include <opencv2/highgui/highgui.hpp>  // imread
#include <opencv2/imgproc/imgproc.hpp>  // cvtColor

#if GF2_USE_PCL
#   include "pcl/common/intersections.h"   // lineWithLineIntersection
#   include "pcl/kdtree/kdtree.h"          // nearestneighboursearch
#   include "pcltools/util.hpp"            // addGaussianNoise, fitLinearPrimitive
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
            //! \brief Main functionality to generate lines from points.
            template <  class       PrimitivePrimitiveAngleFunctorT // concept: energyFunctors.h::PrimitivePrimitiveAngleFunctor
                      , class       PrimitiveContainerT             // concept: std::vector<std::vector<LinePrimitive2>>
                      , class       PointContainerT                 // concept: std::vector<PointPrimitive>
                      , class       PrimitiveT                      = typename PrimitiveContainerT::value_type::value_type  // concept: LinePrimitive2
                      , typename    Scalar                          = typename PrimitiveT::Scalar                           // concept: float
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

    template <  class       _PrimitivePrimitiveAngleFunctorT
              , class       _PrimitiveContainerT
              , class       _PointContainerT
              , class       _PrimitiveT
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
        typedef typename outer_const_iterator::value_type::const_iterator inner_const_iterator;

        if ( out_lines.size() ) std::cerr << "[" << __func__ << "]: " << "warning, out_lines not empty!" << std::endl;
        if ( params.patch_population_limit <= 0 ) { std::cerr << "[" << __func__ << "]: " << "error, popfilter is necessary!!!" << std::endl; return EXIT_FAILURE; }

        // (2) Mix and Filter
        // count patch populations
        GidIntMap populations; // populations[patch_id] = all points with GID==patch_id
        {
            processing::calcPopulations( populations, points );
        }

        // convert local fits
        const _Scalar angle_limit( params.angle_limit / params.angle_limit_div );
        int           nlines = 0;

        int l0,l1,l2,l3;
        l0 = l1 = l2 = l3 = 0; // linear indices cached
        // OUTER0
        for ( outer_const_iterator outer_it0  = in_lines.begin();
                                   outer_it0 != in_lines.end();
                                 ++outer_it0, ++l0 )
        {
            int gid0 = -1;
            l1 = 0; // reset linear counter
            // INNER1
            for ( inner_const_iterator inner_it0  = containers::valueOf<_PrimitiveT>(outer_it0).begin();
                                       inner_it0 != containers::valueOf<_PrimitiveT>(outer_it0).end();
                                     ++inner_it0, ++l1 )
            {

                _PrimitiveT const& prim0 = containers::valueOf<_PrimitiveT>( inner_it0 );

                // sanity checks
                if ( !l1 )
                    gid0 = prim0.getTag( _PrimitiveT::GID ); // store gid of first member
                else if ( prim0.getTag(_PrimitiveT::GID) != gid0 )
                    std::cerr << "[" << __func__ << "]: " << "Not good, prims under one gid don't have same GID..." << std::endl;

                l2 = l0; // reset linear counter
                // OUTER1
                int outer_offs = std::distance( in_lines.begin(), outer_it0 );
                for ( outer_const_iterator outer_it1  = in_lines.begin() + outer_offs;
                                           outer_it1 != in_lines.end();
                                         ++outer_it1, ++l2 )
                {
                    int gid1 = -1;
                    l3 = l1; // reset linear counter
                    int inner_offs = std::distance( containers::valueOf<_PrimitiveT>(outer_it0).begin(), inner_it0 );
                    // INNER1
                    for ( inner_const_iterator inner_it1  = containers::valueOf<_PrimitiveT>(outer_it1).begin() + inner_offs;
                                               inner_it1 != containers::valueOf<_PrimitiveT>(outer_it1).end();
                                             ++inner_it1, ++l3 )
                    {
                        if ( l3 >= outer_it1->size() )
                        {
                            std::cerr << "l3 " << l3 << " >= " << outer_it1->size() << " outer_it1->size()" << std::endl;
                        }

                        // test
                        if ( (std::distance( in_lines.begin(), outer_it0 ) != l0)
                             || (std::distance( containers::valueOf<_PrimitiveT>(outer_it0).begin(), inner_it0 ) != l1)
                             || (std::distance( in_lines.begin(), outer_it1 ) != l2)
                             || (std::distance( containers::valueOf<_PrimitiveT>(outer_it1).begin(), inner_it1) != l3)
                             )
                        {
                            std::cerr << "working on "
                                      << std::distance( in_lines.begin(), outer_it0 ) << "(" << l0 << "), "
                                      << std::distance( containers::valueOf<_PrimitiveT>(outer_it0).begin(), inner_it0 ) << "(" << l1 << "), "
                                      << std::distance( in_lines.begin(), outer_it1 ) << "(" << l2 << "), "
                                      << std::distance( containers::valueOf<_PrimitiveT>(outer_it1).begin(), inner_it1) << "(" << l3 << ")\n";
                            fflush(stderr);
                        }


                        _PrimitiveT const& prim1 = containers::valueOf<_PrimitiveT>( inner_it1 );
                        if ( !l3 )
                            gid1 = prim1.getTag( _PrimitiveT::GID );

                        const bool same_line = equal2D(l0,l1,l2,l3);
                        if ( same_line != ((outer_it0 == outer_it1) && (inner_it0 == inner_it1)) )
                        {
                            std::cerr << "iterator based same_line will NOT work..." << std::endl;
                        }

                        bool add0 = same_line, // add0: new primitive at location of prim0, with direction from prim1. We need to keep a copy, so true if same_line.
                             add1 = false;     // add1: new primitive at location of prim1, with direction from prim0

                        _Scalar closest_angle = _Scalar( 0 );
                        {
                            _Scalar angdiff   = _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>( prim0, prim1, angles, &closest_angle );
                            bool    close_ang = angdiff < angle_limit;
                            add0 |= (close_ang || (populations[gid0] < params.patch_population_limit)); // prim0 needs to be close to prim1 in angle, or small to copy the dir of prim1.
                            add1 |= (close_ang || (populations[gid1] < params.patch_population_limit)); // prim1 needs to be close to prim0 in angle, or small to copy the dir of prim0.
                        }

                        // filter by population
                        {
                            // don't copy direction to patch0 from too small patch1
                            if ( (!same_line) && (populations[gid1] < params.patch_population_limit) )
                            {
                                //std::cout << "skipping " << patch_id0 << "<--" << patch_id1 << "since pop is " << populations[patch_id1] << std::endl;
                                add0 = false;
                            }

                            // don't copy direction to patch1 from too small patch0
                            if ( (!same_line) && (populations[gid0] < params.patch_population_limit) )
                            {
                                //std::cout << "skipping " << patch_id0 << "-->" << patch_id1 << "since pop is " << populations[patch_id0] << std::endl;
                                add1 = false;
                            }

                            // debug check
                            if ( !populations[gid0] )
                                std::cerr << "[" << __func__ << "]: " << "pop[" << gid0 << "] is " << populations[gid0] << std::endl;
                            // debug check
                            if ( !populations[gid1] )
                                std::cerr << "[" << __func__ << "]: " << "pop[" << gid1 << "] is " << populations[gid1] << std::endl;
                        }

        #warning todo: filter already copied direction ids

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
                            cand0.setTag( _PrimitiveT::DIR_GID, prim1.getTag(_PrimitiveT::GID) );
                            // insert
                            int check = containers::add( out_lines, gid0, cand0 ).getTag( _PrimitiveT::GID );
                            ++nlines;

                            // debug
                            if ( check != cand0.getTag(_PrimitiveT::GID) ) std::cerr << "gid tag copy check failed" << std::endl;
                        }

                        if ( !same_line && add1 ) // add other copy only, if it's not the same line (we don't want duplicates)
                        {
                            _PrimitiveT cand1 = _PrimitiveT( prim1.pos(), dir1 );
                            cand1.setTag( _PrimitiveT::GID    , prim1.getTag(_PrimitiveT::GID) );
                            cand1.setTag( _PrimitiveT::DIR_GID, prim0.getTag(_PrimitiveT::GID) );

                            // copy line from pid1 to pid
                            int check = containers::add( out_lines, gid1, cand1 ).getTag( _PrimitiveT::DIR_GID );
                            ++nlines;

                            // debug
                            if ( check != cand1.getTag(_PrimitiveT::DIR_GID) ) std::cerr << "dirgid tag copy check failed" << std::endl;
                        } //...if add1
                    } //...for l3
                } //...for l2
            } //...for l1
        } //...for l0

        // log
        std::cout << "[" << __func__ << "]: " << "finished generating, we now have " << nlines << " candidates" << std::endl;

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::generate()

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

