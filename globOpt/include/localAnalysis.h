#ifndef __GF2_LOCALANALYSIS_H__
#define __GF2_LOCALANALYSIS_H__

#include <vector>
#include <set>
#include "my_types.h"

#if GF2_USE_PCL
    #include "pcl/visualization/pcl_visualizer.h"
    #include "pcl/common/common.h"
#endif

#include "AMUtil2.h"
#include "lineClustering.hpp"

namespace am
{
    class LocalAnalysis
    {
        public:
            template <typename Scalar> inline static std::pair<Scalar,Scalar>
            matching( Eigen::Matrix<Scalar,3,1> const& p1, Eigen::Matrix<Scalar,3,1> const& p2, std::vector<Scalar> const& desired_angles );
#if 0
            template <class LinesT, class PointsT> static inline int
            run( std::vector<PrimitiveClustering<LinesT> > &clusterings, PointsT cloud, float scale );
#endif

            template <class LinesT, class PointsT, class Scalar = typename LinesT::value_type::Scalar > static inline int
            runSimpler( std::vector<PrimitiveClustering<LinesT> >  & clusterings
                        , PointsT                                    cloud
                        , Scalar                              const  scale_arg                = 0.03f
                        , Scalar                              const  similarity_threshold     = 0.005f
                        , int                                 const  max_neighbourhood_size   = 10
                        , std::vector<Scalar>                 const& desired_angles           = {}
                        , std::pair<Scalar,Scalar> (*matching)(Eigen::Matrix<Scalar,3,1> const& p1, Eigen::Matrix<Scalar,3,1> const& p2, std::vector<Scalar> const& desired_angles) = &matching
                        , int                                 const  addNFriends              = 0
                        , Scalar                              const  filter_coeff             = static_cast<Scalar>(1)
                        , LinesT                                   * input_primitives         = NULL
                        , Scalar                              const  trunc_at                 = 1.f
                        );

            template <class LinesT, class PointsT> static inline int
            propose( LinesT      & lines
                     , PointsT                 cloud
                     , std::vector<int> const* indices
                     , int                     K
                     , float                   radius
                     );

            template <class LinesT, class PointsT> static inline int
            display( LinesT                                &lines
                     , PointsT                              cloud
                     , typename LinesT::value_type::Scalar  scale
                     , PrimitiveClustering<LinesT>                  *p_clustering = NULL
                     , const int disp_limit = 10
                     );

            template <class PointsT> static inline int
            plot( PointsT cloud );

            template <class PrimitivesT, typename Scalar = typename PrimitivesT::value_type::Scalar>
            static inline int
            filter( PrimitiveClustering<PrimitivesT>          & out
                    , PrimitiveClustering<PrimitivesT>  const & in
                    , Scalar                            const   working_scale
                    , Scalar                            const   ang_similarity_threshold );

#           if GF2_USE_PCL
#           endif

    };

} // ns am

#ifndef __GF2_LOCALANALYSIS_HPP_INC__
#   define __GF2_LOCALANALYSIS_HPP_INC__
#   include "localAnalysis.hpp"
#endif

// cache compilation
#if 1

#include "primitives/linePrimitive.h"
#include "pcl/point_cloud.h"
namespace am
{
#if 0
    extern template /*<class TLines, class PointsT>*/
    int
    LocalAnalysis::run( std::vector<PrimitiveClustering<std::vector<LinePrimitive> > > &clusterings, boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB> > cloud, float scale );
#endif

    extern template std::pair<float,float>
    LocalAnalysis::matching( Eigen::Matrix<float,3,1> const& p1, Eigen::Matrix<float,3,1> const& p2, std::vector<float> const& desired_angles );
}
#endif

#endif //__GF2_LOCALANALYSIS_H__
