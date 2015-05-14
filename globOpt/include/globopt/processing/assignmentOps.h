#ifndef GO_ASSIGNMENTOPS_H
#define GO_ASSIGNMENTOPS_H

namespace globopt
{

    // Usage: .../Release/bin/toGlobFit --subsample-primitives 0.1 --pop-limit 100 --prims segments.csv --cloud cloud.ply -a points_segments.csv --scale 0.005
    template < class _PrimitiveVectorT
               , class _PrimitiveMapT
               , class _PointContainerT
               , class _PclCloudT>
    int approxUnassignedWPlanes( int argc, char** argv );

}

#endif // GO_ASSIGNMENTOPS_H
