#ifndef GF2_CONTAINERS_HPP
#define GF2_CONTAINERS_HPP

#include <map>
#include <vector>
#include <set>

namespace GF2 {
namespace containers {

    // vector< vector >
    template <class _InnerT> static inline
    _InnerT& add( std::vector<std::vector<_InnerT> > &prims, int gid, _InnerT const& primitive )
    {
        if ( prims.size() <= gid )
            prims.resize( gid + 1 );
        prims[gid].push_back( primitive );

        return prims[gid].back();
    }

    // map< vector >
    template <class _InnerT> static inline
    _InnerT& add( std::map<int,std::vector<_InnerT> > &prims, int gid, _InnerT const& primitive )
    {
        prims[gid].push_back( primitive );
        return prims[gid].back();
    }

    // map< set >
    template <class _InnerT> static inline
    _InnerT add( std::map<int,std::set<_InnerT> > &prims, int gid, _InnerT const& primitive )
    {
        prims[gid].insert( primitive );
        return ( *prims[gid].find(primitive) );
    }

    // vector< vector >::const_iterator
    template <class _PrimitiveT> static inline
    typename std::vector<_PrimitiveT> const&
    valueOf( typename std::vector<std::vector<_PrimitiveT> >::const_iterator const& it )
    {
        return *it;
    }

    // vector< vector >::iterator
    template <class _PrimitiveT> static inline
    typename std::vector<_PrimitiveT> &
    valueOf( typename std::vector<std::vector<_PrimitiveT> >::iterator & it )
    {
        return *it;
    }

    // map< vector >::const_iterator
    template <class _PrimitiveT> static inline
    std::vector<_PrimitiveT> const& valueOf( typename std::map<int,std::vector<_PrimitiveT> >::const_iterator const& it )
    {
        return it->second;
    }

    // map< vector >::iterator
    template <class _PrimitiveT> static inline
    std::vector<_PrimitiveT> & valueOf( typename std::map<int,std::vector<_PrimitiveT> >::iterator & it )
    {
        return it->second;
    }

    // vector::const_iterator
    template <class _PrimitiveT> static inline
    _PrimitiveT const& valueOf( typename std::vector<_PrimitiveT>::const_iterator const& it )
    {
        return *it;
    }

    // vector::iterator
    template <class _PrimitiveT> static inline
    _PrimitiveT & valueOf( typename std::vector<_PrimitiveT>::iterator & it )
    {
        return *it;
    }

} //...ns containers
} //...ns GF2

#endif // GF2_CONTAINERS_HPP
