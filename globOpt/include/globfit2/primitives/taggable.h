#ifndef __GF2_TAGGABLE_H__
#define __GF2_TAGGABLE_H__

#include <string>
#include <map>

namespace GF2
{
    //! \brief      Super-class for the capability of storing tags (ids) in a vector.
    //!
    //!             The types are hard-coded to speed up compilation,
    //!             and because we don't know how many types we will need later.
    class Taggable
    {
        public:
            enum USER_TAGS {
                  USER_ID1 = 10 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID2 = 11 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID3 = 12 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID4 = 13 //!< additional flag to store processing attributes (values only in the generation scope)
                , USER_ID5 = 14 //!< additional flag to store processing attributes (values only in the generation scope)
            }; //...TAGS

            //! \brief              Stores tag for int key. Mostly used by the enum typedefs.
            //! \param[in] key      Key to store \p value at.
            //! \param[in] value    Value to store.
            //! \return             EXIT_SUCCESS
            inline Taggable&
            setTag( int key, int value )
            {
                _tags[key] = value;
                //return EXIT_SUCCESS;
                return *this;
            }

            //! \brief              Returns tag for int key.
            //! \param[in] key      Key to store \p value at.
            //! \return             The value stored for the key. Default -1, if entry is missing.
            inline int
            getTag( int key ) const
            {
                std::map<int,int>::const_iterator it = _tags.find( key );
                if ( it != _tags.end() )
                    return _tags.at( key );
                else
                    return -1;
            }

            //! \brief              Stores tag for string key.
            //! \param[in] key      Key to store \p value at.
            //! \param[in] value    Value to store.
            //! \return             EXIT_SUCCESS
            inline Taggable&
            setTag( std::string key, int value )
            {
                _str_tags[key] = value;
                return *this;
                //return EXIT_SUCCESS;
            }

            //! \brief              Returns tag for string key.
            //! \param[in] key      Key to store \p value at.
            //! \return             The value stored for the key. Default -1, if entry is missing.
            inline int
            getTag( std::string key ) const
            {
                std::map<std::string,int>::const_iterator it = _str_tags.find( key );
                if ( it != _str_tags.end() )
                    return _str_tags.at( key );
                else
                    return -1;
            }

            inline int
            copyTagsFrom( Taggable const& other )
            {
                _tags     = other._tags;
                _str_tags = other._str_tags;

                return EXIT_SUCCESS;
            }

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

        protected:
            std::map<int,int>           _tags;      //!< \brief Stores values for int keys.
            std::map<std::string,int>   _str_tags;  //!< \brief Stores values for string keys.
    }; // ... cls Taggable
} // ... ns GF2

#endif // __GF2_TAGGABLE_H__

