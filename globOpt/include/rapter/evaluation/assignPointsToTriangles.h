#ifndef GO_ASSIGNPOINTSTOTRIANGLES_H
#define GO_ASSIGNPOINTSTOTRIANGLES_H

namespace rapter
{

    // Usage: .../Release/bin/toGlobFit
    template < class _PrimitiveVectorT
             , class _PrimitiveMapT
             , class _PointContainerT
             , class _PclCloudT>
    int assignPointsToTriangles( int argc, char** argv );

} //...ns rapter


#endif // GO_ASSIGNPOINTSTOTRIANGLES_H
