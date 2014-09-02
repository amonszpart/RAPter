#ifndef __GF2_LINEPRIMITIVE2_H__
#define __GF2_LINEPRIMITIVE2_H__

#include "globfit2/primitives/linePrimitive.h"
#include "globfit2/primitives/taggable.h"

namespace GF2
{
    //! \brief      Taggable LinePrimitive. Separate class for historical reasons.
    class LinePrimitive2 : public ::GF2::LinePrimitive, public ::GF2::Taggable
    {
            typedef ::GF2::LinePrimitive ParentT;
        public:
            //! \brief Defines the tags (ids) that this primitive can manage using setTag and getTag functions.
            enum TAGS {
                GID         //!< group id             - which group this primitive is supposed to explain
                , DIR_GID   //!< direction group id   - which group this primitive got it's direction from
            };//...TAGS

#if __cplusplus > 199711L
            //! \brief Inheriting contructor.
            using ::GF2::LinePrimitive::LinePrimitive;
#else
            LinePrimitive2() : ParentT() {}

            //! \brief Constructor that takes raw data in Eigen format as input.
            LinePrimitive2( Eigen::Matrix<Scalar,Dim,1> coeffs ) : ParentT( coeffs ) {}

            //! \brief Constructor that takes raw data in std::vector format as input.
            LinePrimitive2( std::vector<Scalar> const& coeffs ) : ParentT( coeffs ) {}

            LinePrimitive2( Eigen::Matrix<Scalar,3,1> const& p0, Eigen::Matrix<Scalar,3,1> const& dir ) : ParentT( p0, dir ) {}
#endif

    }; //... class LinePrimitive2

} // ... ns GF2

#endif // __GF2_LINEPRIMITIVE2_H__
