#include "rapter/globOpt_types.h"
#include "rapter/primitives/linePrimitive.h"
#include "rapter/primitives/impl/linePrimitive.hpp"

namespace rapter
{
    template
    class Primitive<2,6, rapter::Scalar>;

    LinePrimitive line;

    namespace templ_inst
    {
        typedef Eigen::CwiseUnaryOp<Eigen::internal::scalar_cast_op<double, float>,
                                    Eigen::Map<Eigen::Matrix<double, 3, 1, 0, 3, 1>, 0, Eigen::Stride<0, 0> > const
                                   > DoubleVector3CastAsFloat;

        typedef std::vector<PidT> IndicesContainerT;
    }

    template int
    LinePrimitive::generateFrom( LinePrimitive                              & out
                               , templ_inst::DoubleVector3CastAsFloat  const& normal
                               , float                                 const  distanceFromOrigin );
    template int
    LinePrimitive::generateFrom( LinePrimitive                              & out
                               , LinePrimitive::Position               const& normal
                               , float                                 const  distanceFromOrigin );

    template int
    LinePrimitive::getExtent< PointPrimitiveT
                            , templ_inst::IndicesContainerT
                            , PointContainerT
                            >( LinePrimitive::ExtremaT            & minMax
                             , PointContainerT               const& cloud
                             , double                        const  threshold
                             , templ_inst::IndicesContainerT const* indices_arg
                             , bool                          const  force_axis_aligned ) const;

    template int
    LinePrimitive::draw<PointPrimitiveT, templ_inst::IndicesContainerT, PointContainerT>
                       ( LinePrimitive                    const& line
                       , PointContainerT                  const& cloud
                       , LinePrimitive::Scalar            const  radius
                       , templ_inst::IndicesContainerT    const* indices
                       , pcl::visualization::PCLVisualizer::Ptr  v
                       , std::string                      const  plane_name
                       , double                           const  r
                       , double                           const  g
                       , double                           const  b
                       , int                              const  viewport_id
                       , LinePrimitive::Scalar            const  stretch
                       , int                              const  /*draw_mode*/
                       , float                            const  /*hull_alpha*/
                       );
} //...ns rapter
