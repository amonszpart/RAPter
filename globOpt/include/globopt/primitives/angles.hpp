#ifndef GO_ANGLES_HPP2
#define GO_ANGLES_HPP2

#include "globfit2/simple_types.h" //...__Scalar

namespace globopt
{
    struct AnglesT : public std::vector<GF2::__Scalar>
    {
        typedef GF2::__Scalar       Scalar;
        typedef std::vector<Scalar> ParentT;
        AnglesT() {}
        AnglesT( std::vector<Scalar> const& angles ) : ParentT( angles ) {}
    };

} //...ns globopt

namespace GF2
{
    typedef globopt::AnglesT AnglesT;
}

#endif // GO_ANGLES_HPP2
