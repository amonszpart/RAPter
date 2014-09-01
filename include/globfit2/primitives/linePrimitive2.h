#ifndef __GF2_LINEPRIMITIVE2_H__
#define __GF2_LINEPRIMITIVE2_H__

#include "globfit2/primitives/linePrimitive.h"
#include "globfit2/primitives/taggable.h"

namespace GF2
{
    //! \brief      Taggable LinePrimitive. Separate class for historical reasons.
    class LinePrimitive2 : public ::GF2::LinePrimitive, public ::GF2::Taggable
    {
        public:
            //! \brief Defines the tags (ids) that this primitive can manage using setTag and getTag functions.
            enum TAGS {
                GID         //!< group id             - which group this primitive is supposed to explain
                , DIR_GID   //!< direction group id   - which group this primitive got it's direction from
            };//...TAGS

            //! \brief Inheriting contructor.
            using ::GF2::LinePrimitive::LinePrimitive;

    }; //... class LinePrimitive2

} // ... ns GF2

#endif // __GF2_LINEPRIMITIVE2_H__
