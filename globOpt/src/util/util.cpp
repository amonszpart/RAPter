#include "globfit2/util/util.hpp"

namespace GF2 {

std::string util::timestamp2Str()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time ( &rawtime );
    timeinfo = localtime (&rawtime);

    strftime ( buffer, 80, "_%Y%m%d_%H%M", timeinfo );

    return std::string( buffer );
}

} // namespace GF2
