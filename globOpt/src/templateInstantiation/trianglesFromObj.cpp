#include <vector>
#include "globopt/primitives/impl/triangle.hpp"
#include "globopt/io/impl/trianglesFromObj.hpp"
#include "globfit2/simple_types.h"

namespace globopt
{
    namespace templ_inst
    {
        typedef Triangle<GF2::__Scalar> Triangle;
        typedef std::vector<Triangle>   TrianglesContainer;
    }

    template
    int getTrianglesFromObj( templ_inst::TrianglesContainer & triangles, std::string const& meshPath );
}
