#include "Eigen/Dense"

#include "globopt/primitives/impl/planePrimitive.hpp"
#include "globfit2/globOpt_types.h"

GF2::PlanePrimitive
GF2::PlanePrimitive::fromFileEntry( std::vector<GF2::PlanePrimitive::Scalar> const& entries )
{
    return GF2::PlanePrimitive( /*    pos: */ Eigen::Map<const Eigen::Matrix<Scalar,3,1> >( entries.data()  , 3 ),
                                /* normal: */ Eigen::Map<const Eigen::Matrix<Scalar,3,1> >( entries.data()+3, 3 ) );
}

namespace GF2
{
    std::string PlanePrimitive::toFileEntry() const
    {
        char line[1024];
        sprintf( line, "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,"
                 , pos   ()(0), pos   ()(1), pos   ()(2)
                 , normal()(0), normal()(1), normal()(2) );
        return std::string( line );
    } //...PlanePrimitive::toFileEntry

    PlanePrimitive::Scalar
    PlanePrimitive::getFiniteDistance( PlanePrimitive::ExtentsT const& extrema, PlanePrimitive::Position const& pnt ) const
    {
        return MyPointFinitePlaneDistanceFunctor::eval( extrema, *this, pnt );
    }
} //...ns GF2

namespace GF2
{
    template
    class Primitive<3,6,GF2::Scalar>;

    namespace templ_inst
    {
        typedef Eigen::CwiseUnaryOp<Eigen::internal::scalar_cast_op<double, float>, Eigen::Map<Eigen::Matrix<double, 3, 1, 0, 3, 1>, 0, Eigen::Stride<0, 0> > const> DoubleVector3CastAsFloat;
        typedef std::vector<PidT> IndicesContainerT;
    }

    template int
    PlanePrimitive::generateFrom( PlanePrimitive                             & out
                                , templ_inst::DoubleVector3CastAsFloat  const& normal
                                , float                                 const  distanceFromOrigin );


    template int
    PlanePrimitive::generateFrom( PlanePrimitive                             & out
                                , PlanePrimitive::Position              const& normal
                                , float                                 const  distanceFromOrigin );

    template int
    PlanePrimitive::getExtent< PointPrimitiveT
                             , templ_inst::IndicesContainerT
                             , PointContainerT
                             >( PlanePrimitive::ExtremaT                      & minMax
                              , PointContainerT                          const& cloud
                              , double                                   const  threshold
                              , templ_inst::IndicesContainerT            const* indices_arg
                              , bool                                     const  force_axis_aligned ) const;

    template int
    PlanePrimitive::draw< PointPrimitiveT
                        , PointContainerT
                        , std::vector<PidT>
                        >
                        ( PlanePrimitive                        const& plane
                        , PointContainerT                       const& cloud
                        , float                                 const  radius
                        , std::vector<PidT>                     const* indices
                        , pcl::visualization::PCLVisualizer::Ptr       v
                        , std::string                           const& plane_name
                        , double                                const  r
                        , double                                const  g
                        , double                                const  b
                        , int                                   const  viewport_id  /* = 0*/
                        , float                                 const  stretch      /* = Scalar( 1. ) */
                        , int                                   const  draw_mode    /* = 0*/
                        , float                                 const  alpha        /* = 2.*/
                        );

//    PlanePrimitive::PlanePrimitive( Eigen::Matrix<PlanePrimitive::Scalar, 3, 1> const& pnt,
//                                    Eigen::Matrix<PlanePrimitive::Scalar, 3, 1> const& normal )
//    {
//        _coeffs.template head<3>() = pnt;
//        _coeffs.template segment<3>(3) = normal.normalized();
//    }

    Eigen::Matrix<PlanePrimitive::Scalar,3,1>
    PlanePrimitive::projectPoint( Eigen::Matrix<PlanePrimitive::Scalar,3,1> const& point ) const
    {
        return point - (this->getDistance(point) * this->dir() );
    }

    PlanePrimitive plane;
}
