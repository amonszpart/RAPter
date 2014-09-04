#ifndef GF2_PROC_UTIL_HPP
#define GF2_PROC_UTIL_HPP

#include <map>
#include <set>
#include <vector>
#include "globfit2/util/containers.hpp" // add()

namespace GF2 {         //!< \brief globOpt namespace

    typedef std::map< int, int >    GidIntMap;
    typedef std::set< int >         PidSet;
    typedef std::vector< int >      PidVector;
    typedef std::map< int, PidSet > GidPidSetMap;
    typedef std::map< int, PidVector > GidPidVectorMap;

    namespace processing {  //!< \brief Various 3D processing snippets that don't fit elsewhere.

    //! \brief                  Calculate the number of points that are assigned to a GID (group id).
    //! \tparam                 Concept: std::map< int, int >. Key: gid, value: population (\#points assigned to a group id (GID).
    //! \param[out] populations Contains an int value for each GID that occurred in points. The value is the count, how many points had this GID.
    //! \param[in]  points      The points containing gids retrievable by \code points[pid].getTag( _PointPrimitiveT::GID ) \endcode.
    //! \return                 EXIT_SUCCESS
    template < class _GidIntMap
             , class _PointContainerT >
    inline int calcPopulations( _GidIntMap & populations, _PointContainerT const& points )
    {
        typedef typename _PointContainerT::value_type _PointPrimitiveT;

        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            ++populations[ points[pid].getTag(_PointPrimitiveT::GID) ];
        }

        return EXIT_SUCCESS;
    } //...calcPopulation

    //! \brief                  Calculate the number of points that are assigned to a GID (group id).
    //! \tparam                 Concept: std::map< int, std::set<int> >. Key: GID, value: list of point ids that have that GID. (\#points assigned to a group id (GID).
    //! \param[out] populations Contains an int value for each GID that occurred in points. The value is the count, how many points had this GID.
    //! \param[in]  points      The points containing gids retrievable by \code points[pid].getTag( _PointPrimitiveT::GID ) \endcode.
    //! \return                 EXIT_SUCCESS
    template < class _GidIntSetMap
             , class _PointContainerT >
    inline int getPopulations( _GidIntSetMap & populations, _PointContainerT const& points )
    {
        typedef typename _PointContainerT::value_type _PointPrimitiveT;

        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            const int gid = points[pid].getTag( _PointPrimitiveT::GID );
            containers::add( populations, gid, static_cast<int>(pid) );
        }

        return EXIT_SUCCESS;
    } //...getPopulations

} //...ns processing
} //...ns GF2

#endif // GF2_PROC_UTIL_HPP
