#ifndef GF2_PROC_UTIL_HPP
#define GF2_PROC_UTIL_HPP

#include <map>
#include <set>
#include <vector>
#include "globfit2/util/containers.hpp" // add()

/*!   \brief GlobOpt main namespace. */
namespace GF2 {

    typedef std::map   < int, int >         GidIntMap;
    typedef std::set   < int >              PidSet;
    typedef std::vector< int >              PidVector;
    typedef std::map   < int, PidSet >      GidPidSetMap;
    typedef std::map   < int, PidVector >   GidPidVectorMap;

    /*! \brief Various 3D processing snippets that don't fit elsewhere. */
    namespace processing
    {

        /*! \brief Calculate the number of points that are assigned to a GID (group id).
         *
         *  \tparam                 Concept: std::map< int, int >. Key: gid, value: population (\#points assigned to a group id (GID).
         *  \param[out] populations Contains an int value for each GID that occurred in points. The value is the count, how many points had this GID.
         *  \param[in]  points      The points containing gids retrievable by \code points[pid].getTag( _PointPrimitiveT::GID ) \endcode.
         *  \return                 EXIT_SUCCESS
        */
        template <class _GidIntMap, class _PointContainerT > inline int
        calcPopulations( _GidIntMap & populations, _PointContainerT const& points )
        {
            typedef typename _PointContainerT::value_type _PointPrimitiveT;

            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                ++populations[ points[pid].getTag(_PointPrimitiveT::GID) ];
            }

            return EXIT_SUCCESS;
        } //...calcPopulation

        /*! \brief Calculate the number of points that are assigned to a GID (group id).
        *   \tparam                 Concept: std::map< int, std::set<int> >. Key: GID, value: list of point ids that have that GID. (\#points assigned to a group id (GID).
        *   \param[out] populations Contains an int value for each GID that occurred in points. The value is the count, how many points had this GID.
        *   \param[in]  points      The points containing gids retrievable by \code points[pid].getTag( _PointPrimitiveT::GID ) \endcode.
        *   \return                 EXIT_SUCCESS
        */
        template <class _GidIntSetMap, class _PointContainerT > inline int
        getPopulations( _GidIntSetMap & populations, _PointContainerT const& points )
        {
            typedef typename _PointContainerT::value_type _PointPrimitiveT;

            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                const int gid = points[pid].getTag( _PointPrimitiveT::GID );
                containers::add( populations, gid, static_cast<int>(pid) );
            }

            return EXIT_SUCCESS;
        } //...getPopulations

        /*! \brief Adds angles from generator to in/out vector \p angles.
         *
         *  \tparam _AnglesContainerT Concept: std::vector<_Scalar>.
         *  \tparam _Scalar           Floating point precision type.
         *  \param[in,out] angles     Container (vector) to append to.
         *  \param[in] angle_gen      Generator element. 0, \p angle_gen, 2 * \p angle_gen, ..., M_PI will be appended.
         *  \param[in] verbose        Enable logging.
         */
        template <class _AnglesContainerT, class _Scalar> inline int
        appendAnglesFromGenerator( _AnglesContainerT &angles, _Scalar &angle_gen, char verbose = true )
        {
            std::set<_Scalar> angles_set;
            //angles_set.insert( )
            if ( std::find(angles.begin(), angles.end(), _Scalar(0)) != angles.end() )
            {
                // log
                if ( verbose )
                {
                    std::cout << "[" << __func__ << "]: " << "adding 0 to ";
                    for(size_t vi=0;vi!=angles.size();++vi)
                        std::cout << angles[vi] << ((vi==angles.size()-1) ? "" : ", ");
                    std::cout << "\n";
                }

                angles.push_back( _Scalar(0) );
            }

            // generate
            for ( _Scalar angle = angle_gen; angle < M_PI; angle+= angle_gen )
                angles.push_back( angle );

            if ( std::find(angles.begin(), angles.end(), _Scalar(M_PI)) != angles.end() )
            {
                // log
                if ( verbose )
                {
                    std::cout << "[" << __func__ << "]: " << "adding " << M_PI << " to ";
                    for ( size_t vi=0;vi!=angles.size();++vi)
                        std::cout << angles[vi] << ((vi==angles.size()-1) ? "" : ", ");
                    std::cout << "\n";
                }

                angles.push_back( _Scalar(M_PI) );
            }

            // print
            if ( verbose )
            {
                std::cout << "Desired angles: {";
                for ( size_t vi=0;vi!=angles.size();++vi)
                    std::cout << angles[vi] << ((vi==angles.size()-1) ? "" : ", ");
                std::cout << "}\n";
            }

            return EXIT_SUCCESS;
        } //...appendAnglesFromGenerator()

    } //...ns processing
} //...ns GF2

#endif // GF2_PROC_UTIL_HPP
