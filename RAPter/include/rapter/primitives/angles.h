#ifndef RAPTER_ANGLES_HPP
#define RAPTER_ANGLES_HPP

#include "rapter/simpleTypes.h" //...__Scalar

namespace rapter
{
    struct AnglesT : public std::vector<rapter::__Scalar>
    {
        typedef rapter::__Scalar       Scalar;
        typedef std::vector<Scalar> ParentT;
        AnglesT() {}
        AnglesT( std::vector<Scalar> const& angles ) : ParentT( angles ) {}
    };

} //...ns rapter

#endif // RAPTER_ANGLES_HPP
