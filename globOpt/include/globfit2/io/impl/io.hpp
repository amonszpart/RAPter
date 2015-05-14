#include "globfit2/io/io.h"
#include <fstream>

#ifdef WITH_QCQPCPP
#   include "qcqpcpp/io/io.h"
#endif

namespace GF2
{
    namespace io
    {

    } // ... ns io
} // ... ns GF2

namespace globopt
{
    namespace io
    {
        template <class _PrimitiveContainerT, class _PointContainerT>
        inline std::string writePrimAssocCloud( _PrimitiveContainerT const& prims, _PointContainerT const& points, std::string stem, std::string dir )
        {
            typedef typename _PointContainerT::PrimitiveT PointPrimitiveT;

            char outPrimsPath[ 512 ], outAssocPath[512], outCloudPath[512];
            sprintf( outPrimsPath, "%s/%s.csv", dir.c_str(), stem.c_str() );
            std::cout << " writing " << prims.size() << " primitives to " << outPrimsPath << std::endl;
            GF2::io::savePrimitives<typename _PrimitiveContainerT::PrimitiveT,typename _PrimitiveContainerT::mapped_type::const_iterator >( prims, std::string(outPrimsPath) );

            sprintf( outAssocPath, "%s/points_%s.csv", dir.c_str(), stem.c_str() );
            std::cout << " writing " << points.size() << " assignments to " << outAssocPath << std::endl;
            GF2::io::writeAssociations<PointPrimitiveT>( points, outAssocPath );

            sprintf( outCloudPath, "%s/%s.cloud.ply", dir.c_str(), stem.c_str() );
            std::cout << " writing " << points.size() << " points to " << outCloudPath << std::endl;
            GF2::io::writePoints<PointPrimitiveT>( points, outCloudPath );

            std::stringstream ss;
            ss << "../show.py -p " << outPrimsPath << " -a " << outAssocPath << " --cloud " << outCloudPath;

            return ss.str();
        }
    }
}
