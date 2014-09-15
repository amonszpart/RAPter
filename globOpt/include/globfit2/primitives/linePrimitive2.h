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
                GID        = 0  //!< group id             - which group this primitive is supposed to explain
                , DIR_GID  = 1  //!< direction group id   - which group this primitive got it's direction from
                , CHOSEN   = 2  //!< an additional flag to store, if this is part of a solution.
                , USER_ID1 = 10 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID2 = 11 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID3 = 12 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID4 = 13 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID5 = 14 //!< additional flag to store processing attributes (values only in the generation scope)
            };//...TAGS


            typedef ParentT::Scalar Scalar;

#if 0 // __cplusplus > 199711L
            //! \brief Inheriting contructor.
            using ::GF2::LinePrimitive::LinePrimitive;
            //LinePrimitive2() : ParentT() {}
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
