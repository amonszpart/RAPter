#ifndef GF2_PARSE_H
#define GF2_PARSE_H

#if GF2_USE_PCL

#include "pcl/console/parse.h"
#include <vector>
#include <string>
#include "boost/algorithm/string.hpp"


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

    inline int
    parse_argument (int argc, char** argv, const char* str, long int &val)
    {
        int index = pcl::console::find_argument (argc, argv, str) + 1;

        if (index > 0 && index < argc )
          val = atoi (argv[index]);

        return (index - 1);
    }

    inline int parse_x_arguments (int argc, char** argv, const char* str, std::vector<long>& v)
    {
      for (int i = 1; i < argc; ++i)
      {
        // Search for the string
        if ((strcmp (argv[i], str) == 0) && (++i < argc))
        {
          // look for ',' as a separator
          std::vector<std::string> values;
          boost::split (values, argv[i], boost::is_any_of (","), boost::token_compress_on);

          v.resize (values.size ());
          for (size_t j = 0; j < v.size (); ++j)
            v[j] = atoi (values.at (j).c_str ());

          return (i - 1);
        }
      }
      return (-1);
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
