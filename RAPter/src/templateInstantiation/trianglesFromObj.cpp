#include <vector>
#include "rapter/primitives/impl/triangle.hpp"
#include "rapter/io/impl/trianglesFromObj.hpp"
#include "rapter/simpleTypes.h"

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
