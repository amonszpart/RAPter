#ifndef GF2_DISKUTIL_HPP
#define GF2_DISKUTIL_HPP

#include <string>
#include "boost/filesystem.hpp"

namespace GF2 {
namespace util {

//! \brief  _saveSafe Saves backups of the path to be overwritten. If .0.bak exists, try .1.bak, etc.
//! \param            Path, that will be overwritten
//! \return           Path saved
inline std::string
saveBackup( std::string path )
{
    std::string bak_path = path;
    if ( boost::filesystem::exists(path) )
    {
        int i = -1;
        do
        {
            ++i;
            std::stringstream ss;
            ss << path << "." << i << ".bak";
            bak_path = ss.str();

        } while ( boost::filesystem::exists(bak_path) );

        std::cout << "saving backup of " << path << " to " << bak_path << std::endl;
        boost::filesystem::copy( path, bak_path );
    }

    return bak_path;
} // ..._saveSafe()

} //...namespace util
} //...namespace GF2

#endif // GF2_DISKUTIL_HPP
