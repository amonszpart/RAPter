#ifndef GF2_CONTAINERS_HPP
#define GF2_CONTAINERS_HPP

#include <map>
#include <vector>
#include <set>
#include "globfit2/simple_types.h"
#include "globfit2/util/exception.h"

namespace GF2 {
namespace containers {

    // vector< vector >
    template <class _InnerT> static inline
    _InnerT& add( std::vector<std::vector<_InnerT> > &prims, GidT gid, _InnerT const& primitive )
    {
        if ( prims.size() <= gid )
            prims.resize( gid + 1 );
        prims[gid].push_back( primitive );

        return prims[gid].back();
    }

    // map< vector >
    template <class _InnerT> static inline
    _InnerT& add( std::map<GidT,std::vector<_InnerT> > &prims, GidT gid, _InnerT const& primitive )
    {
        prims[gid].push_back( primitive );
        return prims[gid].back();
    }

    // map< set >
    template <class _InnerT> static inline
    _InnerT add( std::map<GidT,std::set<_InnerT> > &prims, GidT gid, _InnerT const& primitive )
    {
        prims[gid].insert( primitive );
        return ( *prims[gid].find(primitive) );
    }

    // vector
    template <class _InnerT> static inline
    _InnerT& add( std::vector<_InnerT> &container, _InnerT const& value )
    {
        container.push_back( value );

        return container.back();
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
    std::vector<_PrimitiveT> const& valueOf( typename std::map<GidT,std::vector<_PrimitiveT> >::const_iterator const& it )
    {
        return it->second;
    }

    // map< vector >::iterator
    template <class _PrimitiveT> static inline
    std::vector<_PrimitiveT> & valueOf( typename std::map<GidT,std::vector<_PrimitiveT> >::iterator & it )
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


    /*!
     * \tparam _StdIteratorType      Either map<GidT,vector<PrimitiveT> >::iterator or map<>::const_iterator
     * \tparam _StdInnerIteratorType Either vector<PrimitiveT>::iterator or vector<>::const_iterator
     */
    template <class _ParentT, class _PrimitiveT, class _StdIteratorType, class _InnerStdIteratorType>
    struct PrimitiveContainerIterator
    {
        protected:
            inline _InnerStdIteratorType startInnerAt( _StdIteratorType &outer_it )
            {
                _lid1 = 0;
                return outer_it->second.begin();
            }

        public:
            GENERATE_CLASS_EXCEPTION("PrimitiveContainerIterator")

            inline PrimitiveContainerIterator( _ParentT& container )
                : _container( container )
                , _outer_it ( _container.begin() )
                , _inner_it ( startInnerAt( _outer_it ) )
                , _lid0( 0 )
                , _lid1( 0 )
                , _uid( 0 )
            {}

            inline void step()
            {
                if ( _outer_it == _container.end() )
                {
                    //std::cerr << "[PrimitiveContainer::Iterator] Iterators stepped beyond container's end!" << std::endl;
                    throw PrimitiveContainerIterator::Exception("Iterators stepped beyond container's end!");
                }

                bool step_outer = false;
                // we have valid _outer_it, so step inner
                if ( _inner_it != _outer_it->second.end() )
                {
                    ++_inner_it;
                    ++_lid1;
                    if ( _inner_it == _outer_it->second.end() )
                        step_outer = true;
                }
                else // step outer
                {
                    step_outer = true;
                }

                if ( step_outer )
                {
                    ++_outer_it;
                    ++_lid0;

                    // if outer valid, restart inner
                    if ( _outer_it != _container.end() )
                        _inner_it = startInnerAt( _outer_it );
                }

                ++_uid;
            } //...step()

            /*! \brief Returns gid of current primitive.
             *  \return _outer_it->first, which should be equal to _inner_it.getTag( _PrimitiveT::TAGS::GID ).
             */
            inline GidT getGid() const
            {
                if ( _inner_it->getTag( _PrimitiveT::TAGS::GID ) != _outer_it->first )
                    throw PrimitiveContainerIterator::Exception("_outer_it->first != prim.GID ");

                return _outer_it->first;
            }

            inline DidT getDid() const
            {
                return _inner_it->getTag( _PrimitiveT::TAGS::DIR_GID );
            }

            /*! \brief Returns linear index of gid in the map of gids.
             */
            inline LidT getLid0() const
            {
                return _lid0;
            }

            /*! \brief Returns linear index of primitive in the vector of primitives contained under the current gid key in the map.
             */
            inline LidT getLid1() const
            {
                return _lid1;
            }

            inline bool hasNext() const
            {
                if ( _outer_it != _container.end() && (_inner_it == _outer_it->second.end()) )
                    throw PrimitiveContainerIterator::Exception("[hasNext()] illegal state, _inner_it is at inner_end.");

                return _outer_it != _container.end();
            }

            inline bool operator<( PrimitiveContainerIterator const& other )
            {
                //"[PrimitiveContainer::Iterator] HACK, assumes GIDs come in ascending order! TODO: change to iterator comparison"
                if   ( getGid() < other.getGid() )  return true;
                else                                return (getLid0() < other.getLid0());
            }

            //inline _PrimitiveT* operator->() const
            inline typename std::iterator_traits<_InnerStdIteratorType>::pointer operator->() const
            {
                if ( _outer_it == _container.end() || (_inner_it == _outer_it->second.end()) )
                    throw PrimitiveContainerIterator::Exception("[Iterator] operator-> called of invalid iterator!");

                return _inner_it.operator->();
            }

            //inline _PrimitiveT& operator*() const
            inline typename std::iterator_traits<_InnerStdIteratorType>::reference operator*() const
            {
                if ( _outer_it == _container.end() || (_inner_it == _outer_it->second.end()) )
                    throw PrimitiveContainerIterator::Exception("[Iterator] operator* called of invalid iterator!");

                return *_inner_it;
            }

            inline UidT getUniqueId() const
            {
                return _uid;
            }

        protected:
            _ParentT                            &_container;
            _StdIteratorType                    _outer_it;
            _InnerStdIteratorType               _inner_it;
            LidT                                _lid0, _lid1;
            UidT                                _uid;

    }; //...struct Iterator

    template <class _PrimitiveT>
    class PrimitiveContainer : public std::map< GidT, std::vector<_PrimitiveT> >
    {
        public:
            GENERATE_CLASS_EXCEPTION("PrimitiveContainer")

            typedef _PrimitiveT                              PrimitiveT;
            typedef std::vector<_PrimitiveT>                 InnerContainerT;
            typedef std::map< GidT, InnerContainerT >        ParentT;
            typedef typename ParentT::iterator               ParentIteratorT;
            typedef typename InnerContainerT::iterator       InnerContainerIteratorT;
            typedef typename ParentT::const_iterator         ParentConstIteratorT;
            typedef typename InnerContainerT::const_iterator InnerContainerConstIteratorT;

            typedef PrimitiveContainerIterator<ParentT,_PrimitiveT,ParentIteratorT, InnerContainerIteratorT>          Iterator;
            typedef PrimitiveContainerIterator<ParentT,_PrimitiveT,ParentConstIteratorT,InnerContainerConstIteratorT> ConstIterator;
    }; //...struct PrimitiveContainer


} //...ns containers
} //...ns GF2

#endif // GF2_CONTAINERS_HPP
