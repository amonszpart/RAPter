#ifndef RAPTER_TAGGABLE_HPP__
#define RAPTER_TAGGABLE_HPP__

#include <string>
#include <map>
#include "rapter/primitives/taggable.h"
#include "rapter/simpleTypes.h"

namespace rapter
{
    /*! \brief              Stores long tag for int key. Mostly used by the enum typedefs.
      * \param[in] key      Key to store \p value at.
      * \param[in] value    Value to store.
      * \return             EXIT_SUCCESS
      */
    template <typename _Scalar>
    inline Taggable<_Scalar>&
    Taggable<_Scalar>::setTag( int key, long value )
    {
        _longValuedTags[key] = value;
        return *this;
    }

    /*! \brief              Stores tag for int key.
     *  \param[in] key      Key to store \p value at.
     *  \param[in] value    Value to store.
     *  \return             EXIT_SUCCESS
     */
    template <typename _Scalar>
    inline Taggable<_Scalar>&
    Taggable<_Scalar>::setTag( int key, int value )
    {
        _longValuedTags[key] = value;
        return *this;
    }

    /*! \brief              Stores size_t valued tag for int key.
     *  \param[in] key      Key to store \p value at.
     *  \param[in] value    Value to store.
     *  \return             EXIT_SUCCESS
     */
    template <typename _Scalar>
    inline Taggable<_Scalar>&
    Taggable<_Scalar>::setTag( int key, std::size_t value )
    {
        _longValuedTags[key] = value;
        return *this;
    }

    /*! \brief              Stores tag for int key.
     *  \param[in] key      Key to store \p value at.
     *  \param[in] value    Value to store.
     *  \return             EXIT_SUCCESS
     */
    template <typename _Scalar>
    inline Taggable<_Scalar>&
    Taggable<_Scalar>::setTag( char key, char value )
    {
        _charValuedTags[key] = value;
        return *this;
    }

    /*! \brief              Stores tag for string key.
     *  \param[in] key      Key to store \p value at.
     *  \param[in] value    Value to store.
     *  \return             EXIT_SUCCESS
     */
    template <typename _Scalar>
    inline Taggable<_Scalar>&
    Taggable<_Scalar>::setTag( _Scalar key, _Scalar value )
    {
        _scalarValuedTags[key] = value;
        return *this;
    }

    // ______________________________________________________________________

    /*! \brief              Returns tag for int key.
     *  \param[in] key      Key to store \p value at.
     *  \return             The value stored for the key. Default -1, if entry is missing.
     */
    template <typename _Scalar>
    inline rapter::GidT
    Taggable<_Scalar>::getTag( rapter::GidT key ) const
    {
        typename std::map<rapter::GidT, rapter::GidT>::const_iterator it = _longValuedTags.find( key );
        if ( it != _longValuedTags.end() )
            return it->second; // _tags.at( key ); changed by Aron on 3/1/2015

        return TAG_UNSET;
    }

    template <typename _Scalar>
    inline int
    Taggable<_Scalar>::getTag( int key ) const
    {
        return getTag( static_cast<rapter::GidT>(key) );
    }

    template <typename _Scalar>
    inline char
    Taggable<_Scalar>::getTag( char key ) const
    {
        typename std::map<char,char>::const_iterator it = _charValuedTags.find( key );
        if ( it != _charValuedTags.end() )
            return it->second;

        return TAG_UNSET;
    }

    /*! \brief              Returns tag for int key.
     *  \param[in] key      Key to store \p value at.
     *  \return             The value stored for the key. Default -1, if entry is missing.
     */
    template <typename _Scalar>
    inline _Scalar
    Taggable<_Scalar>::getTag( _Scalar key ) const
    {
        typename std::map<_Scalar,_Scalar>::const_iterator it = _scalarValuedTags.find( key );
        if ( it != _scalarValuedTags.end() )
            return it->second; // _tags.at( key ); changed by Aron on 3/1/2015

        return TAG_UNSET;
    }

    template <typename _Scalar>
    inline int
    Taggable<_Scalar>::copyTagsFrom( Taggable<_Scalar> const& other )
    {
        _longValuedTags    = other._longValuedTags;
        _scalarValuedTags  = other._scalarValuedTags;
        _charValuedTags    = other._charValuedTags;
        //_intValuedTags     = other._intValuedTags;
        //_str_tags          = other._str_tags;

        return EXIT_SUCCESS;
    }
} // ... ns rapter

#endif // __GO_TAGGABLE_HPP__

