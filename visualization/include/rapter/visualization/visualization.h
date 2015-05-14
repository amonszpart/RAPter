#ifndef GF2_VISUALIZATION_H
#define GF2_VISUALIZATION_H

#include <vector>
#include "pcl/visualization/pcl_visualizer.h"

#include "globfit2/primitives/pointPrimitive.h"

namespace GF2 {
namespace vis {

    //! \brief Return type for visualizer. Convenience reasons (might need to be changed to std::shared_ptr.
    typedef pcl::visualization::PCLVisualizer::Ptr MyVisPtr;

    //typedef GF2::LinePrimitive2                             PrimitiveT;
    //typedef std::vector< std::vector<PrimitiveT> >          PrimitiveContainerT;
    typedef GF2::PointPrimitive                             PointPrimitiveT;
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
} //...ns GF2

#ifndef GF2_INC_VISUALIZATION_HPP
#   define GF2_INC_VISUALIZATION_HPP
#   include "globfit2/visualization/impl/visualization.hpp"
#endif // GF2_INC_VISUALIZATION_HPP

#endif // GF2_VISUALIZATION_H
