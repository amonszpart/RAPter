#ifndef __GO_TAGGABLE_H__
#define __GO_TAGGABLE_H__

#include <string>
#include <map>
#include <globfit2/simple_types.h>

namespace globopt
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
            Taggable&
            setTag( int key, long value );

            /*! \brief              Stores tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            Taggable&
            setTag( int key, int value );

            /*! \brief              Stores size_t valued tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            Taggable&
            setTag( int key, std::size_t value );

            /*! \brief              Stores tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            Taggable&
            setTag( char key, char value );

            /*! \brief              Stores tag for string key.
             *  \param[in] key      Key to store \p value at.
             *  \param[in] value    Value to store.
             *  \return             EXIT_SUCCESS
             */
            Taggable&
            setTag( _Scalar key, _Scalar value );

            // ______________________________________________________________________

            /*! \brief              Returns tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \return             The value stored for the key. Default -1, if entry is missing.
             */
            GF2::GidT
            getTag( GF2::GidT key ) const;

            int
            getTag( int key ) const;

            char
            getTag( char key ) const;

            /*! \brief              Returns tag for int key.
             *  \param[in] key      Key to store \p value at.
             *  \return             The value stored for the key. Default -1, if entry is missing.
             */
            _Scalar
            getTag( _Scalar key ) const;

            int
            copyTagsFrom( Taggable const& other );


        protected:
            std::map<GF2::GidT,GF2::GidT>         _longValuedTags;     //!< \brief Stores long,int and char values for int keys.
            std::map<char,char>         _charValuedTags;
            std::map<_Scalar,_Scalar>   _scalarValuedTags;   //!< \brief Stores values for int keys with value float/double.
            //std::map<std::string,int>   _str_tags;           //!< \brief Stores values for string keys.
    }; // ... cls Taggable
} // ... ns globopt

#endif // __GO_TAGGABLE_H__


