#ifndef GF2_GLOBFIT_H
#define GF2_GLOBFIT_H

namespace globopt
{
    template <class _PclCloudT, class _PrimitiveMapT, class _PrimitiveVectorT, class _PointContainerT>
    int toGlobFit( int argc, char **argv );

    template <class _PclCloudT, class _PrimitiveMapT, class _PrimitiveVectorT, class _PointContainerT>
    int fromGlobFit( int argc, char **argv );
} //...ns GF2

#endif // GF2_GLOBFIT_H
