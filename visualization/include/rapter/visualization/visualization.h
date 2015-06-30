#ifndef RAPTER_VISUALIZATION_H
#define RAPTER_VISUALIZATION_H

#include <vector>
#include "pcl/visualization/pcl_visualizer.h"

#include "rapter/primitives/pointPrimitive.h"

namespace rapter {
namespace vis {

    //! \brief Return type for visualizer. Convenience reasons (might need to be changed to std::shared_ptr.
    typedef pcl::visualization::PCLVisualizer::Ptr MyVisPtr;

    //typedef rapter::LinePrimitive2                             PrimitiveT;
    //typedef std::vector< std::vector<PrimitiveT> >          PrimitiveContainerT;
    typedef rapter::PointPrimitive                             PointPrimitiveT;
    typedef std::vector<PointPrimitiveT>                    PointContainerT;
    typedef typename PointPrimitiveT::Scalar                Scalar;

    //! \brief Preimplemented, C++03 version for lines that reads info from disk.
    template <class PrimitiveT> int
    showCli( int argc, char** argv );

    //! \brief Preimplemented, C++03 version for lines.
//        vis::MyVisPtr
//        showLines( vis::lines::PrimitiveContainerT  const& lines
//                 , vis::lines::PointContainerT      const& points
//                 , vis::lines::Scalar               const  scale    = lines::Scalar(0.05)
//                 , Eigen::Matrix<lines::Scalar,3,1> const& colour   = (Eigen::Matrix<lines::Scalar,3,1>() << 0,0,1).finished()
//                 , bool                             const  spin     = true
//                 , std::vector<lines::Scalar>       const* angles   = NULL
//                 , bool                             const  show_ids = false
//                 , char                             const  use_tags = true  );


} //...ns vis
} //...ns rapter

#include "rapter/visualization/impl/visualization.hpp"

#endif // RAPTER_VISUALIZATION_H
