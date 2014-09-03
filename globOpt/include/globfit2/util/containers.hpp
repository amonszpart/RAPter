#ifndef GF2_CONTAINERS_HPP
#define GF2_CONTAINERS_HPP

#include <map>
#include <vector>

namespace GF2 {
namespace containers {
    template <class _PrimitiveT> static inline
    _PrimitiveT& add( std::vector<std::vector<_PrimitiveT> > &prims, int gid, _PrimitiveT const& primitive )
    {
        if ( prims.size() <= gid )
            prims.resize( gid + 1 );
        prims[gid].push_back( primitive );

        return prims[gid].back();
    }

    template <class _PrimitiveT> static inline
    _PrimitiveT& add( std::map<int,std::vector<_PrimitiveT> > &prims, int gid, _PrimitiveT const& primitive )
    {
        prims[gid].push_back( primitive );
        return prims[gid].back();
    }

    template <class _PrimitiveT> static inline
    typename std::vector<_PrimitiveT> const&
    valueOf( typename std::vector<std::vector<_PrimitiveT> >::const_iterator const& it )
    {
        return *it;
    }

    template <class _PrimitiveT> static inline
    std::vector<_PrimitiveT> const& valueOf( typename std::map<int,std::vector<_PrimitiveT> >::const_iterator const& it )
    {
        return it->second;
    }

    template <class _PrimitiveT> static inline
    _PrimitiveT const& valueOf( typename std::vector<_PrimitiveT>::const_iterator const& it )
    {
        return *it;
    }

} //...ns containers
} //...ns GF2

#endif // GF2_CONTAINERS_HPP
