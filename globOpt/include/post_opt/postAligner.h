#ifndef __GF2_POSTALIGNER_H__
#define __GF2_POSTALIGNER_H__

#include <vector>
#include <Eigen/Dense>

namespace GF2
{
    class PostAligner
    {
        public:
            template < class      PrimitivesT
                       , class    PointsPtrT
                       , typename Scalar        = typename PrimitivesT::value_type::Scalar
                       , int        Dim
                     > static inline int
            run( PrimitivesT             & prims
                 , PointsPtrT              points
                 , std::vector<int> const& points_prims
                 , std::vector<int> const& prims_clusters
                 , Eigen::Matrix<Scalar,3,1> (*getDirection)(typename PrimitivesT::value_type const& prim) );

            // get line inliers
            template< class PrimitivesT, class PointsPtrT
                      , typename Scalar = typename PrimitivesT::value_type::Scalar
                      , int Dim         = PrimitivesT::value_type::Dim >
            static inline int
            points2Prims( std::vector<int>                      & points_lines
                          , PointsPtrT                            cloud
                          , PrimitivesT                    const& prims
                          , float (*distance)(Eigen::Matrix<Scalar,3,1> const& pnt, Eigen::Matrix<Scalar,Dim,1> const& prim) = NULL );

    };

} // ... ns GF2

#ifndef __GF2_POSTALIGNER_HPP_INC__
#   define __GF2_POSTALIGNER_HPP_INC__
#   include "post_opt/postAligner.hpp"
#endif


#endif // __GF2_POSTALIGNER_H__
