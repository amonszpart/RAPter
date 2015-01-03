#ifndef __GF2_TAGGABLE_H__
#define __GF2_TAGGABLE_H__

#include <string>
#include <map>
#include <globfit2/simple_types.h>

namespace GF2
{
    /*! \brief     Super-class for the capability of storing tags (ids) in a vector.
     *
     *             The types are hard-coded to speed up compilation,
     *             and because we don't know how many types we will need later.
     */
    template <typename _Scalar>
    class Taggable
    {
        public:
            typedef _Scalar Scalar;
            //typedef std::map<int,_Scalar> IntScalarMapT;

            static const int TAG_UNSET = -1;

            enum USER_TAGS {
                  USER_ID1 = 10 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID2 = 11 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID3 = 12 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID4 = 13 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID5 = 14 //!< additional flag to store processing attributes (values only in the generation scope)
            }; //...TAGS

            /*! \brief              Stores long tag for int key. Mostly used by the enum typedefs.
              * \param[in] key      Key to store \p value at.
              * \param[in] value    Value to store.
              * \return             EXIT_SUCCESS
              */
            inline Taggable&
            setTag( int key, long value )
            {
                _longValuedTags[key] = value;
                return *this;
            }

            /*! \brief              Stores tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            inline Taggable&
            setTag( int key, int value )
            {
                _longValuedTags[key] = value;
                return *this;
            }

            /*! \brief              Stores size_t valued tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            inline Taggable&
            setTag( int key, std::size_t value )
            {
                _longValuedTags[key] = value;
                return *this;
            }

            /*! \brief              Stores tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            inline Taggable&
            setTag( char key, char value )
            {
                _charValuedTags[key] = value;
                return *this;
            }

            /*! \brief              Stores tag for string key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            inline Taggable&
            setTag( _Scalar key, _Scalar value )
            {
                _scalarValuedTags[key] = value;
                return *this;
            }

            // ______________________________________________________________________

            /*! \brief              Returns tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \return             The value stored for the key. Default -1, if entry is missing.
             */
            inline GidT
            getTag( GidT key ) const
            {
                typename std::map<GidT,GidT>::const_iterator it = _longValuedTags.find( key );
                if ( it != _longValuedTags.end() )
                    return it->second; // _tags.at( key ); changed by Aron on 3/1/2015

                return TAG_UNSET;
            }

            inline int
            getTag( int key ) const
            {
                return getTag( static_cast<GidT>(key) );
            }

            inline char
            getTag( char key ) const
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
            inline _Scalar
            getTag( _Scalar key ) const
            {
                typename std::map<_Scalar,_Scalar>::const_iterator it = _scalarValuedTags.find( key );
                if ( it != _scalarValuedTags.end() )
                    return it->second; // _tags.at( key ); changed by Aron on 3/1/2015

                return TAG_UNSET;
            }

            inline int
            copyTagsFrom( Taggable const& other )
            {
                _longValuedTags    = other._longValuedTags;
                _scalarValuedTags  = other._scalarValuedTags;
                _charValuedTags    = other._charValuedTags;
                //_intValuedTags     = other._intValuedTags;
                //_str_tags          = other._str_tags;

                return EXIT_SUCCESS;
            }

#if 0
            /*! \brief              Returns tag for string key.
             *  \param[in] key      Key to store \p value at.
             *  \return             The value stored for the key. Default -1, if entry is missing.
             */
            inline int
            getTag( std::string key ) const
            {
                std::map<std::string,int>::const_iterator it = _str_tags.find( key );
                if ( it != _str_tags.end() )
                    return _str_tags.at( key );
                else
                    return TAG_UNSET;
            }

            /*! \brief              Stores tag for string key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            inline Taggable&
            setTag( std::string key, int value )
            {
                _str_tags[key] = value;
                return *this;
            }

            // not safe anymore
           /*! \brief Hack to return all int keyed tags. Used in \ref io::savePrimitives(). TODO: change to Taggable.serialize, deserialize.
             */
            inline std::vector<int>
            getIntTags() const
            {
                std::vector<int> tags;
                for ( std::map<int,int>::const_iterator it = _tags.begin(); it != _tags.end(); ++it )
                    tags.push_back( it->second );
                return tags;
            }
#endif

        protected:
            std::map<GidT,GidT>         _longValuedTags;     //!< \brief Stores long,int and char values for int keys.
            std::map<char,char>         _charValuedTags;
            std::map<_Scalar,_Scalar>   _scalarValuedTags;   //!< \brief Stores values for int keys with value float/double.
            //std::map<std::string,int>   _str_tags;           //!< \brief Stores values for string keys.
    }; // ... cls Taggable
} // ... ns GF2

#endif // __GF2_TAGGABLE_H__

