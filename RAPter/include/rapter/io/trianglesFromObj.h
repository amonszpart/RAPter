#ifndef GO_TRIANGLESFROMOBJ_H
#define GO_TRIANGLESFROMOBJ_H

namespace rapter
{
    template <class _TrianglesContainer, typename Scalar>
    int getTrianglesFromObj( _TrianglesContainer & triangles, std::string const& meshPath, Scalar minPlaneEdge = 1. );
}


#endif // GO_TRIANGLESFROMOBJ_H
