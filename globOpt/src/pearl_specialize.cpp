#if USE_PEARL && GF2_USE_PCL

#include "pearl/pearl.h"

#include  <vector>
#include "pcl/point_cloud.h"
#include "pcl/point_types.h"
#include "primitives/linePrimitive.h"
#include "primitives/planePrimitive.h"

template int
am::Pearl::run(
        std::vector<int>          & labels
        , std::vector<am::LinePrimitive>      & lines
        , boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB> >            const& cloud
        , std::vector<int>   const* indices
        , am::Pearl::Params                    params
        , std::vector<std::vector<int> > *label_history
        , std::vector<std::vector<am::LinePrimitive> > *line_history
        );

template int
am::Pearl::run(
        std::vector<int>          & labels
        , std::vector<am::PlanePrimitive>      & lines
        , boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB> >            const& cloud
        , std::vector<int>   const* indices
        , am::Pearl::Params                    params
        , std::vector<std::vector<int> > *label_history
        , std::vector<std::vector<am::PlanePrimitive> > *line_history
        );

#endif
