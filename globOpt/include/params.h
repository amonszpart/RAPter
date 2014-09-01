#ifndef __GF2_PARAMS_H__
#define __GF2_PARAMS_H__

#include <map>

namespace GF2
{
    template <typename T>
    class Params
    {
            typedef std::map<std::string,T> MapType;

        public:
            inline T get( std::string key ) const
            {
                typename MapType::const_iterator it = _map.find( key );
                if ( it != _map.end() )
                {
                    return it->second;
                }

                return T(-1);
            }

            inline bool has( std::string key ) const
            {
                return _map.find( key ) != _map.end();
            }

            inline int set( std::string key, T value )
            {
                _map[ key ] = value;
                return 0;
            }

        protected:
             MapType _map;
    };
}

#endif // __GF2_PARAMS_H__
