#ifndef GO_SUBSAMPLEPRIMITIVES_HPP
#define GO_SUBSAMPLEPRIMITIVES_HPP

#include "globfit2/simple_types.h"
#include "globfit2/io/inputParser.hpp"
#include "globfit2/processing/util.hpp"

namespace globopt
{
    template < class _PrimitiveVectorT
             , class _PrimitiveMapT
             , class _PointContainerT
             , class _PclCloudT>
    int subsamplePrimitives( int argc, char** argv )
    {
        typedef typename _PrimitiveMapT::mapped_type InnerContainerT;
        typedef typename _PclCloudT::Ptr             PclPtrT;
        typedef typename InnerContainerT::value_type PrimitiveT;
        typedef typename PrimitiveT::Scalar Scalar;
        typedef typename _PointContainerT::PrimitiveT PointPrimitiveT;

        _PointContainerT    points;
        PclPtrT             pclCloud;
        _PrimitiveVectorT   primitivesVector;
        _PrimitiveMapT      primitives;
        struct MyParams { Scalar scale; } params;

        GF2::parseInput<InnerContainerT,_PclCloudT>( points, pclCloud, primitivesVector, primitives, params, argc, argv, /* readAssoc: */ true );

        Scalar ratio( 1. );
        if ( GF2::console::parse_argument(argc,argv,"--subsample-primitives",ratio) < 0 )
        {
            std::cerr << "[" << __func__ << "]: " << "you need to provide subsample ratio after --subsample-primitives\n";
            return EXIT_FAILURE;
        }


        GidPidVectorMap populations;
        processing::getPopulations( populations, points );
        for ( int )
    } //...subsamplePrimitives()

}

#endif // GO_SUBSAMPLEPRIMITIVES_HPP
