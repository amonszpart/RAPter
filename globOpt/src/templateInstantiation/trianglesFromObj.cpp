#include <vector>
#include "globopt/primitives/impl/triangle.hpp"
#include "globopt/io/impl/trianglesFromObj.hpp"
#include "globfit2/simple_types.h"

namespace rapter
{
    namespace templ_inst
    {
        typedef Triangle<rapter::__Scalar> Triangle;
        typedef std::vector<Triangle>   TrianglesContainer;
    }

    template
    int getTrianglesFromObj( templ_inst::TrianglesContainer & triangles, std::string const& meshPath, rapter::__Scalar minPlaneEdge );

}
