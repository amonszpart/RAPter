#ifndef __GF2_PEARL_H__
#define __GF2_PEARL_H__

#include <vector>
#include <Eigen/Dense>
#include "pcl/point_cloud.h"
#include <limits>

namespace am
{
    class LinePrimitive;

    class Pearl
    {
        public:
            struct Params
            {
                    Params()
                        : scale                 ( 0.02f )
                        , beta                  ( int_mult * 100 )
                        , lambdas               ( 4,1 )
                        , gammasqr              ( 250.f )
                        , int_mult              ( 1e3f )
                        , data_max              ( std::numeric_limits<int>::max() / 2 ) // TODO: change back to 1e10, if breaks
                        , max_neighbourhood_size( 15 )
                        , max_pearl_iterations  ( 3 )
                        , pottsweight           ( 1.f )
                        , debug                 ( true )
                    {
                        lambdas << 0.f, 0.f, 1.f, 0.f;
                    }

                    float   scale;                     //!< \brief "max_neighbourhood_radius", analogous to globOpt.scale.
                    float   beta;                      //!< \brief Complexity cost.
                    Eigen::Matrix<float,-1,1> lambdas; //!< \brief Pairwise cost: lambdas(2) * exp( -1.f * distsqr / gammasqr ). Other parts of lambdas not used.
                    float   gammasqr;                  //!< \brief Pairwise cost: lambdas(2) * exp( -1.f * distsqr / gammasqr );

                    float   int_mult;                  //!< \brief Scales up distances to show up in "int" domain.
                    int     data_max;                  //!< \brief Data term maximum (to prevent int overflow of global energy)

                    int     max_neighbourhood_size;    //!< \brief Sync with \ref CandidateGeneratorParams::nn_K.
                    int     max_pearl_iterations;      //!< \brief How many iterations of pearl to do.

                    float   pottsweight;               //!< \brief Unused.
                    bool    debug;
            };

            template <class _PclCloudT, class _PointsContainerT, class _PrimitiveT>
            static inline int
            run( std::vector<int>                       & labels
               , std::vector<_PrimitiveT>               & lines
               , _PclCloudT                        const& cloud
               , _PointsContainerT                 const& points
               , std::vector<int>                  const* indices
               , Params                            const& params = Params()
               , std::vector<std::vector<int        > > *label_history = NULL
               , std::vector<std::vector<_PrimitiveT> > *line_history = NULL
               , int                               const nPropose = 0 // how many primitives to propose, if none provided
               );

            template <typename PointT, class TLine>
            static inline int
            propose( std::vector<TLine> &lines
                     , boost::shared_ptr<pcl::PointCloud<PointT> > cloud
                     , std::vector<int> const* indicies
                     , Params             const& params
                     , int              const nPropose = 0
                    );

            template <class _PointsContainerT, typename _PclPointT, class _PrimitiveT, typename Scalar = float>
            static inline int
            expand( std::vector<int>          & labels
                    , boost::shared_ptr<pcl::PointCloud<_PclPointT> > cloud
                    , _PointsContainerT   const& points
                    , std::vector<int>   const* indices
                    , std::vector<_PrimitiveT> const& lines
                    , Scalar             const  gammasqr
                    , Scalar             const  beta
                    , Params             const& params
                    );

            template <class PointsT, class TLine>
            static inline int
            refit( std::vector<TLine>        &lines
                   , std::vector<int> const& labels
                   , PointsT          const& cloud
//                   , std::vector<int> const* indices
                   , Params             const& params
                   );

            inline int
            getActiveLabelsCount( std::vector<int>         const& labels
                                  , std::vector<unsigned>       * active_labels_arg = NULL );
    };
} // nsam

#ifndef __GF2_INC_PEARL_HPP__
#   define __GF2_INC_PEARL_HPP__
#   include "pearl.hpp"
#endif

#endif // __GF2_PEARL_H__
