#ifndef GF2_MERGING_H
#define GF2_MERGING_H

#include <iostream>

namespace GF2 {

class Merging
{
    public:

        template < class    _PrimitiveContainerT
                 , class    _PointContainerT
                 , typename _Scalar          = float
                 , class    _PointPrimitiveT = typename _PointContainerT::value_type
                 , class    _PrimitiveT      = typename _PrimitiveContainerT::value_type::value_type
                 >
        static inline int mergeCli( int argc, char** argv );

};

} //... namespace GF2


#ifndef GF2_INC_MERGING_HPP
#   define GF2_INC_MERGING_HPP
#   include "globfit2/optimization/impl/merging.hpp"
#endif

#endif // GF2_MERGING_H
