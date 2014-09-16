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
            }; //...TAGS

            typedef ParentT::Scalar Scalar;

            LinePrimitive2() : ParentT() {}

            //! \brief Constructor that takes raw data in Eigen format as input.
            LinePrimitive2( Eigen::Matrix<Scalar,Dim,1> coeffs ) : ParentT( coeffs ) {}

            //! \brief Constructor that takes raw data in std::vector format as input.
            LinePrimitive2( std::vector<Scalar> const& coeffs ) : ParentT( coeffs ) {}

            LinePrimitive2( Eigen::Matrix<Scalar,3,1> const& p0, Eigen::Matrix<Scalar,3,1> const& dir ) : ParentT( p0, dir ) {}

            LinePrimitive2( Eigen::Matrix<Scalar,3,1> const& centroid, Eigen::Matrix<Scalar,3,1> const& eigen_values, Eigen::Matrix<Scalar, 3, 3> const& eigen_vectors )
            : ParentT( centroid, eigen_values, eigen_vectors ) {}

            // _______________________IO_______________________
            /*! \brief Used in \ref io::readPrimitives to determine how many floats to parse from one entry.
             */
            static inline int getFileEntryLength() { return 6; }

            /*! \brief Output <x0,n>, location and normal for the plane to a string that does *not* have an endline at the end.
             */
            inline std::string toFileEntry() const
            {
                Eigen::Matrix<Scalar,3,1> pos = this->pos();
                Eigen::Matrix<Scalar,3,1> nrm = this->normal();
                char line[1024];
                sprintf( line, "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,"
                         , pos(0), pos(1), pos(2)
                         , nrm(0), nrm(1), nrm(2) );
                return std::string( line );
            }

            /*! \brief              Used in \ref io::readPrimitives to create the primitive from the read floats.
             *  \param[in] entries  Contains <x0,n> to create the plane from.
             */
            static inline LinePrimitive2 fromFileEntry( std::vector<Scalar> const& entries )
            {
                return LinePrimitive2( Eigen::Map<const Eigen::Matrix<Scalar,3,1> >( entries.data()  , 3 ),
                                       Eigen::Map<const Eigen::Matrix<Scalar,3,1> >( entries.data()+3, 3 ).cross(Eigen::Matrix<Scalar,3,1>::UnitZ()) );
            } //...fromFileEntry()


    }; //... class LinePrimitive2

} // ... ns GF2

#endif // __GF2_LINEPRIMITIVE2_H__
