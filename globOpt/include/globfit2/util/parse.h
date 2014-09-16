#ifndef GF2_PARSE_H
#define GF2_PARSE_H

#if GF2_USE_PCL

#include "pcl/console/parse.h"


namespace GF2 {

// Here follow a bunch of relay functions that will make it easier to detach the project from PCL.
// Todo: move to util.hpp or parse.hpp

namespace console
{
    /*! \brief Find switch.
     */
    inline bool find_switch( int argc, char** argv, const char* argument_name )
    {
        return pcl::console::find_switch( argc, argv, argument_name );
    }

    template <typename T>
    inline int parse_argument (int argc, char** argv, const char* str, T &val)
    {
        return pcl::console::parse_argument( argc, argv, str, val );
    }

//    inline int parse_argument (int argc, char** argv, const char* str, bool &val)
//    {
//        pcl::console::parse_argument( argc, argv, val );
//    }

//    inline int parse_argument (int argc, char** argv, const char* str, float &val)
//    {
//        pcl::console::parse_argument( argc, argv, val );
//    }

//    inline int parse_argument (int argc, char** argv, const char* str, double &val)
//    {
//        pcl::console::parse_argument( argc, argv, val );
//    }
} //...ns console
} //...ns GF2

#endif //...GF2_USE_PCL

#endif // GF2_PARSE_H
