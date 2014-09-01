#ifndef GLOBFIT2_H
#define GLOBFIT2_H

#include <string>

#include "pcl/point_cloud.h"
#include "pcl/ModelCoefficients.h"
#include "pcl/visualization/pcl_visualizer.h"
#include "my_types.h"
#include "kmeans/primitive.h"
#include "primitives/linePrimitive.h"

#include "opencv2/core.hpp"

#define M_PIf (static_cast<float>(M_PI))

namespace am
{

    template <class TLines> class PrimitiveClustering;

    class GlobFit2
    {
        public:
            //typedef std::vector<LinePrimitive> TLines;

            static int
            image_2_2DCloud( pcl::PointCloud<MyPoint>::Ptr      & cloud
                             , cv::Mat                     const& img
                             , int                                N_samples = 100
                             , float                       const  Z         = -1.f
                             , float                       const  scale     = 1.f );
            static int
            image_2_2DCloud( pcl::PointCloud<MyPoint>::Ptr &cloud
                             , std::string                 img_path
                             , int                         N_samples = 100
                             , const float                 Z         = -1.f
                             , const float                 scale     = 1.f );

            template <class TLines>
            static int
            optimize( MaskType                                  & min_config
                      , TLines                             const& lines
                      , PrimitiveClustering<TLines>        const* clustering
                      , pcl::PointCloud<MyPoint>::Ptr             cloud
                      , Eigen::Matrix<float,-1,1>          const& lambdas
                      , float                              const  scale          = .03f
                      , std::vector<float>                 const& desired_angles = {0, M_PI_2, M_PI}
                      , int                                const  max_iterations = 16
                      , int                                const  max_step_count = 700
                      , double                             const  threshold      = 0.01
                      , float                              const  trunc_pw_at_angle = .25f // angle in radians
                      , std::vector<int>                   const* indices        = NULL
                      , int                                const  fixedK         = -1
                      , std::vector<int>                        * labels         = NULL
                      , float                                   * min_e_arg      = NULL
                      , std::string                        const* p_out_dir      = NULL );

            template <class PrimitivesT, class PointsPtrT, typename Scalar = typename PrimitivesT::value_type::Scalar> static inline int
            assignPoints( std::vector<int>          * points_lines
                          //, GF2::SelectionType const& selection
                          , MaskType           const* mask
                          , PrimitivesT        const& primitives
                          , PointsPtrT                cloud
                          , std::vector< std::vector<int> > *primitives_points = NULL
                          , Scalar             const  dist_threshold           = -1.f // to assign points to primitives
                          , std::vector<Scalar>     * p_distances_sqr          = NULL );

            template <class PrimitivesT, class PointsPtrT, typename Scalar = typename PrimitivesT::value_type::Scalar> static inline int
            split( std::vector<std::vector<Scalar> > &lines_ends
                   , MaskType    const& mask
                   , PrimitivesT const& prims
                   , PointsPtrT         cloud
                   , Scalar      const  scale
                   , std::vector<int>  * points_lines_arg );

            /**
             * @brief visualize Shows input and output of optimize.
             * @param lines
             * @param optimized_mask
             * @param real_mask
             * @param cloud
             * @param threshold
             * @param indices
             * @param fixedK
             * @return
             */
            template <class PrimitivesT> static int
            visualize( PrimitivesT                          const& lines
                       , std::vector<int>                   const& optimized_mask
                       , pcl::PointCloud<MyPoint>::Ptr             cloud
                       , const float                               scale     = 0.01
                       , std::vector<float>                 const& desired_angles = { 0, M_PI_2, M_PI }
                       , MaskType                           const* real_mask = NULL
                       , pcl::PointIndices::Ptr                    indices   = pcl::PointIndices::Ptr()
                       , std::vector<int>                   const* labels    = NULL
                       , bool                               const  spin      = true
                       , bool                               const  show_angles = false
                       , float                              const  gap_threshold = 0.2f );

            static int
            scoreParallelity( std::vector<int> mask
                              , std::vector<LinePrimitive> const& lines );

            template <class PrimitivesT, typename Scalar> static inline Scalar
            queryEnergy( Eigen::Matrix<Scalar,-1,1>                & energies
                         , MaskType                           const& mask
                         , PrimitivesT                        const& lines
                         , pcl::PointCloud<MyPoint>::ConstPtr        cloud
                         , Eigen::Matrix<Scalar,-1,1>         const& lambdas
                         , Scalar                             const  scale
                         , std::vector<Scalar>                const& desired_angles
                         , double                             const  threshold
                         , float                              const  trunc_at_pw_angle
                         , std::vector<int>                   const* indices = NULL
                         );

    }; // class globfit2

} // ns am

#ifndef __GF2_GLOBFIT2_HPP_INC__
#   define __GF2_GLOBFIT2_HPP_INC__
#   include "globFit2.hpp"
#endif

#endif // GLOBFIT2_H
