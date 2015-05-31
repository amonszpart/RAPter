#ifndef RAPTER_ASSIGNPOINTSTOTRIANGLES_H
#define RAPTER_ASSIGNPOINTSTOTRIANGLES_H

namespace rapter
{

    // Usage: .../Release/bin/toGlobFit
    template < class _PrimitiveVectorT
             , class _PrimitiveMapT
             , class _PointContainerT
             , class _PclCloudT>
    int assignPointsToTriangles( int argc, char** argv );

} //...ns rapter


#endif // RAPTER_ASSIGNPOINTSTOTRIANGLES_H
