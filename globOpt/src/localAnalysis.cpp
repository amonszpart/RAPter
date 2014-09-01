#include "localAnalysis.h"
#include "primitives/linePrimitive.h"
#include "primitives/planePrimitive.h"

namespace am
{
#if 0
    //template <class TLines, class PointsT>
    template int
    LocalAnalysis::run( std::vector<PrimitiveClustering<std::vector<LinePrimitive> > > &clusterings, boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB> > cloud, float scale );
#endif

    template std::pair<float,float>
    LocalAnalysis::matching( Eigen::Matrix<float,3,1> const& p1, Eigen::Matrix<float,3,1> const& p2, std::vector<float> const& desired_angles );

    //template <class LinesT, class PointsT> inline
    template int
    LocalAnalysis::runSimpler(
            std::vector<PrimitiveClustering<std::vector<LinePrimitive> > >      & clusterings
            , boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB> > cloud
            , float                              const  scale
            , float                              const  similarity_threshold
            , int                                const  max_neighbourhood_size
            , std::vector<float>                 const& desired_angles           = {}
            , std::pair<float,float> (*matching)(Eigen::Matrix<float,3,1> const& p1, Eigen::Matrix<float,3,1> const& p2, std::vector<float> const& desired_angles) = &(matching<float>)
            , int                                const  addNFriends              = 0
            , float                              const  filter                   = static_cast<float>(1)
            , std::vector<LinePrimitive>              * input_primitives         = NULL
            , float                              const  trunc_at                = 1.f
            );

    template int
    LocalAnalysis::runSimpler(
            std::vector<PrimitiveClustering<std::vector<PlanePrimitive> > >      & clusterings
            , boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB> > cloud
            , float                              const  scale
            , float                              const  similarity_threshold
            , int                                const  max_neighbourhood_size
            , std::vector<float>                 const& desired_angles           = {}
            , std::pair<float,float> (*matching)(Eigen::Matrix<float,3,1> const& p1, Eigen::Matrix<float,3,1> const& p2, std::vector<float> const& desired_angles) = &(matching<float>)
            , int                                const  addNFriends              = 0
            , float                              const  filter                   = static_cast<float>(1)
            , std::vector<PlanePrimitive>              * input_primitives        = NULL
            , float                              const  trunc_at                = 1.f
            );
} // nsam
