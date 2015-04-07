#ifndef GF2_GLOBFIT_HPP
#define GF2_GLOBFIT_HPP
#include "globfit2/comparison/globFit.h"

#include <iostream>
#include "globfit2/io/inputParser.hpp"
#include "globfit2/util/util.hpp"       // timeStamp2Str
#include "globfit2/util/containers.hpp"

namespace globopt
{
#if 0
    bool loadGlobfit( const std::string& filename )
    {
        std::string line;
        std::ifstream fin(filename.c_str());

        if ( !fin.good() )
        {
            std::cout << "[" << __func__ << "]: " << "!fin.good() -> return false;" << std::endl;
            return false;
        }

        size_t numPoint = 0;
        while ( fin )
        {
            getline(fin, line);
            if (line.empty() || line[0]=='#' || line[0] == 13)
            {
                continue;
            }

            std::stringstream sin(line);
            sin >> numPoint;
            std::cout << "[" << __func__ << "]: " << "numPoints: " << numPoint << std::endl;
            break;
        }

        while (fin && numPoint != 0)
        {
            getline(fin, line);

            if (line.empty() || line[0]=='#' || line[0] == 13)
            {
                continue;
            }

            --numPoint;

            RichPoint* p = new RichPoint();
            _vecPointSet.push_back(p);

            double x, y, z, nx, ny, nz, confidence;

            std::stringstream sin(line);
            sin >> x >> y >> z >> nx >> ny >> nz >> confidence;

            p->point        = Point(x, y, z);
            p->normal       = Vector(nx, ny, nz)/sqrt(nx*nx+ny*ny+nz*nz);
            p->confidence   = confidence;
        }

        size_t numPrimitive = 0;
        while (fin)
        {
            getline(fin, line);
            if ( line.empty() || line[0]=='#' || line[0] == 13)
            {
                std::cout << "[" << __func__ << "]: " << "skipping line " << line << std::endl;
                continue;
            }

            std::cout << "[" << __func__ << "]: " << "parsing line \"" << line << "\" == " << (int)line[0] << std::endl;
            std::stringstream sin(line);
            sin >> numPrimitive;
            std::cout << "[" << __func__ << "]: " << "numPrimitive: " << numPrimitive << std::endl;
            break;
        }

        while (fin && numPrimitive != 0){
            getline(fin, line);

            if (line.empty() || line[0]=='#' || line[0] == 13)
            {
                continue;
            }

            std::string indication = line.substr(0, line.find_first_of(" \t"));
            Primitive* pPrimitive = NULL;
            if (indication == "plane") {
                pPrimitive = new Plane(_vecPointSet);
            } else if (indication == "cylinder") {
                pPrimitive = new Cylinder(_vecPointSet);
            } else if (indication == "cone") {
                pPrimitive = new Cone(_vecPointSet);
            } else if (indication == "sphere") {
                pPrimitive = new Sphere(_vecPointSet);
            }

            if (pPrimitive == NULL) {
                std::cerr << "Error: bad file format!" << std::endl;
                continue;
            }

            fin.seekg(-((int)(line.size()+1)), std::ios::cur);
            if (pPrimitive->load(fin))
            {
                pPrimitive->setIdx(_vecPrimitive.size());
                _vecPrimitive.push_back(pPrimitive);
                --numPrimitive;
                //std::cout << "[" << __func__ << "]: " << "read primitive " << std::endl;
            } else {
                delete pPrimitive;
                std::cerr << "Error: bad file format!" << std::endl;
            }
        }

        return true;
    }
#endif

    template <class _PclCloudT, class _PrimitiveMapT, class _PrimitiveVectorT, class _PointContainerT>
    int fromGlobFit( int argc, char **argv )
    {
        typedef typename _PclCloudT::Ptr                            PclCloudPtrT;
        typedef typename _PrimitiveMapT::InnerContainerT            InnerPrimitiveContainerT;
        typedef typename _PrimitiveMapT::PrimitiveT                 PrimitiveT;
        typedef typename PrimitiveT::Scalar                         Scalar;
        typedef typename GF2::containers::PrimitiveContainer<PrimitiveT> IterablePrimitiveContainer;
        typedef typename _PointContainerT::PrimitiveT               PointPrimitiveT;
        typedef GF2::PidT PidT;

        std::cout << "hello fromGlobfit" << std::endl;
        std::string primitivesPath("");
        if ( GF2::console::parse_argument( argc, argv, "--from", primitivesPath ) < 0 )
        {
            std::cerr << "[" << __func__ << "]: " << "need --from primitivesPath to exist!" << std::endl;
            return EXIT_FAILURE;
        }
        std::string inputSegmentsPath; // the original ransac input
        if (    (GF2::console::parse_argument( argc, argv, "-p"     , inputSegmentsPath) < 0)
                && (GF2::console::parse_argument( argc, argv, "--prims", inputSegmentsPath) < 0)
                && (!boost::filesystem::exists(inputSegmentsPath)) )
        {
            std::cerr << "[" << __func__ << "]: " << "-p or --prims is compulsory" << std::endl;
            return EXIT_FAILURE;
        }

        Scalar scale;
        if (    (GF2::console::parse_argument( argc, argv, "-s"     , scale) < 0)
             && (GF2::console::parse_argument( argc, argv, "--scale", scale) < 0)
                )
        {
                std::cerr << "need scale, -s or --scale" << std::endl;
                return EXIT_FAILURE;
        }

        std::string outName( "primitives" );
        GF2::console::parse_argument( argc, argv, "-o", outName);
        std::cout << "saving to " << outName << ".globfit.csv, change with -o if needed" << std::endl;

        _PointContainerT        points;
        //PclCloudPtrT            pclCloud;
        _PrimitiveVectorT       /*primitivesVector,*/ segmentsVector;
        _PrimitiveMapT          primitives, segments;
        //bool                    valid_input = true;
        int                     err = EXIT_SUCCESS;

        std::cout << "[" << __func__ << "]: " << "reading primitives from " << inputSegmentsPath << "...";
        err = GF2::io::readPrimitives<PrimitiveT, InnerPrimitiveContainerT>( segmentsVector, inputSegmentsPath, &segments );
        if ( EXIT_SUCCESS != err )
        {
            std::cerr << "could not read input segments..." << std::endl;
            return err;
        }

        std::string   line;
        std::ifstream inStream( primitivesPath.c_str() );

        if ( !inStream.good() )
        {
            std::cout << "[" << __func__ << "]: " << "!fin.good() -> return false;" << std::endl;
            return false;
        }

        PidT numPoints = 0;
        while ( inStream )
        {
            getline( inStream, line );
            if ( line.empty() || line[0]=='#' || line[0] == 13 )
                continue;

            std::stringstream sin( line );
            sin >> numPoints;
            std::cout << "[" << __func__ << "]: " << "numPoints: " << numPoints << std::endl;
            break;
        }

        points.resize( numPoints );
        while ( inStream && numPoints )
        {
            getline( inStream, line );

            if ( line.empty() || line[0]=='#' || line[0] == 13 )
                continue;

            --numPoints;

            double x, y, z, nx, ny, nz, confidence;

            std::stringstream sin( line );
            sin >> x >> y >> z >> nx >> ny >> nz >> confidence;

            points.at( numPoints ).setTag( PointPrimitiveT::TAGS::PID, numPoints );
//            std::cout << "setting pid " << numPoints
//                      << ": " << points.at( numPoints ).getTag( PointPrimitiveT::TAGS::PID )
//                      << "\tfrom line " << line << std::endl;
            points.at( numPoints ).template coeffs().template head<3>   ( ) = Eigen::Matrix<Scalar,3,1>(x,y,z);
            points.at( numPoints ).template coeffs().template segment<3>(3) = Eigen::Matrix<Scalar,3,1>(nx,ny,nz).normalized();
            //p->point        = Point(x, y, z);
            //p->normal       = Vector(nx, ny, nz)/sqrt(nx*nx+ny*ny+nz*nz);
            //p->confidence   = confidence;
        }

        size_t numPrimitives = 0;
        while ( inStream )
        {
            getline( inStream, line );
            if ( line.empty() || line[0]=='#' || line[0] == 13)
                continue;

            std::stringstream sin( line );
            sin >> numPrimitives;
            //std::cout << "[" << __func__ << "]: " << "numPrimitive: " << numPrimitives << std::endl;
            break;
        }

        auto segmentsIt = segments.begin();
        PidT planeId = -1;
        while ( inStream && numPrimitives )
        {
            getline( inStream, line );

            if ( line.empty() || line[0]=='#' || line[0] == 13 )
                continue;

            std::string primitiveType = line.substr(0, line.find_first_of(" \t"));

            //Primitive* pPrimitive = NULL;
            if ( primitiveType == "points" )
            {
                std::stringstream sin( line.substr(line.find_first_of(" \t")+1) );
                PidT pid;
                while ( sin >> pid )
                {
                    points.at( pid ).setTag( PointPrimitiveT::TAGS::GID, primitives[planeId].at(0).getTag( PrimitiveT::TAGS::GID ) );
                }
                --numPrimitives;
                ++segmentsIt;
                continue;
            }
            else if (primitiveType == "plane")
            {
                std::stringstream sin( line.substr(line.find_first_of(" \t")+1) );
                Eigen::Matrix<Scalar,3,1> normal;
                Scalar distance;
                sin >> normal(0) >> normal(1) >> normal(2) >> distance;
                //std::cout << "normal: " << normal.transpose() << ", d:" << distance
                //          << ", from line: " << line
                //          << std::endl;
                planeId = numPrimitives-1;
                primitives[ planeId ].push_back( PrimitiveT() );
                PrimitiveT& plane = primitives[planeId].back();
                PrimitiveT::generateFrom( plane, normal, distance );
                //std::cout << "plane: " << plane.toString() << std::endl;
                GF2::GidT gid = segmentsIt->second[0].getTag(PrimitiveT::TAGS::GID);
                GF2::DidT did = segmentsIt->second[0].getTag(PrimitiveT::TAGS::DIR_GID );
                plane.setTag( PrimitiveT::TAGS::GID, gid );
                plane.setTag( PrimitiveT::TAGS::DIR_GID, did );
                plane.setTag( PrimitiveT::TAGS::STATUS, segmentsIt->second[0].getTag(PrimitiveT::TAGS::STATUS)  );
                //std::cout << "equiv input gid was " << segmentsIt->second[0].getTag(PrimitiveT::TAGS::GID) << std::endl;
            }
//            else if (primitiveType == "cylinder")
//            {
//                //pPrimitive = new Cylinder(_vecPointSet);
//                //std::cerr << "[" << __func__ << "]: " << "can't read cylinder" << std::endl;
//            } else if (primitiveType == "cone")
//            {
//                //std::cerr << "[" << __func__ << "]: " << "can't read cone" << std::endl;
//                //pPrimitive = new Cone(_vecPointSet);
//            } else if (primitiveType == "sphere")
//            {
//                //std::cerr << "[" << __func__ << "]: " << "can't read sphere" << std::endl;
//                //pPrimitive = new Sphere(_vecPointSet);
//            }
        } //...while primitives

        std::stringstream ssPrims; ssPrims << outName << ".globfit.csv";
        GF2::io::savePrimitives<PrimitiveT,typename InnerPrimitiveContainerT::const_iterator>( primitives, ssPrims.str() );
        std::stringstream ssAssoc; ssAssoc << "points_" << outName << ".globfit.csv";
        GF2::io::writeAssociations<PointPrimitiveT>( points, ssAssoc.str() );
        std::cout << "results written to " << ssPrims.str() << " and " << ssAssoc.str() << "\n";
        std::cout << "../show.py -s " << scale << " -p " << ssPrims.str() << " -a " << ssAssoc.str() << " --save-poly" << std::endl;

        return EXIT_SUCCESS;
    } //...toGlobFit()

    template <typename _Scalar>
    struct GlobFitParams
    {
        _Scalar scale;
    };

    template <class _PclCloudT, class _PrimitiveMapT, class _PrimitiveVectorT, class _PointContainerT>
    int toGlobFit( int argc, char **argv )
    {
        typedef typename _PclCloudT::Ptr                            PclCloudPtrT;
        typedef typename _PrimitiveMapT::InnerContainerT            InnerContainerT;
        typedef typename _PrimitiveMapT::PrimitiveT                 PrimitiveT;
        typedef typename PrimitiveT::Scalar                         Scalar;
        typedef typename GF2::containers::PrimitiveContainer<PrimitiveT> IterablePrimitiveContainer;
        typedef typename _PointContainerT::PrimitiveT               PointPrimitiveT;
        typedef GF2::PidT PidT;


        std::cout << "hello globfit" << std::endl;

        _PointContainerT        points;
        PclCloudPtrT            pclCloud;
        _PrimitiveVectorT       primitivesVector;
        _PrimitiveMapT          primitives;
        bool                    valid_input = true;
        GlobFitParams<Scalar>   params;

        // parse
        {
            valid_input =
                    (EXIT_SUCCESS ==
                     GF2::parseInput<InnerContainerT,_PclCloudT>( points, pclCloud, primitivesVector, primitives, params, argc, argv ) );
            if ( !valid_input )
            {
                std::cerr << "[" << __func__ << "]: " << "could not parse input..." << std::endl;
                return EXIT_FAILURE;
            }
        }
        std::cout << "have " << primitivesVector.size() << " primitives" << std::endl;
        std::string timeStamp = GF2::util::timestamp2Str();
        std::map< GF2::PidT, std::vector<GF2::PidT> > primitivesPoints;

        std::string     outPath( std::string("segments") /*+ timeStamp*/ + ".globfit" );
        std::ofstream   f      ( outPath );
        f << "# Number of Points\n";
        f << points.size() << std::endl;
        f << "# Here come the " << points.size() << " Points\n";
        f << "# point_x point_y point_z normal_x normal_y normal_z confidence\n";
        GF2::PidT pointId( 0 );
        for ( typename _PointContainerT::const_iterator it = points.begin(); it != points.end(); ++it, ++pointId )
        {
            f << it->template pos()(0) << " " << it->template pos()(1) << " " << it->template pos()(2) << " "
              << it->template dir()(0) << " " << it->template dir()(1) << " " << it->template dir()(2)
              << " 1.0"
              << std::endl;
            primitivesPoints[ it->getTag(PointPrimitiveT::TAGS::GID) ].push_back( pointId );
        }
        f << "# End of Points\n";

        PidT primitiveCount( 0 );
        for ( typename IterablePrimitiveContainer::Iterator it(primitives); it.hasNext(); it.step(), ++primitiveCount ) ;

        f << "\n# Number of Primitives\n";
        f << primitiveCount << std::endl;
        f << "# Here come the " << primitiveCount << " Primitives\n";

        GF2::PidT           primitiveId( 0 );
        std::vector<Scalar> coeffs;
        for ( typename IterablePrimitiveContainer::Iterator it(primitives); it.hasNext(); it.step(), ++primitiveId )
        {
            it->to4Coeffs( coeffs );
            f << "# Primitive " << primitiveId << "\n";
            f << "# plane normal_x normal_y normal_z d\n";
            f << "plane " << coeffs[0] << " " << coeffs[1] << " " << coeffs[2] << " " << coeffs[3] << std::endl;
            f << "# points idx_1 idx_2 idx_3 ...\n";
            f << "points ";
            auto indicesIt = primitivesPoints.find( it->getTag(PrimitiveT::TAGS::GID) );
            if ( indicesIt != primitivesPoints.end() )
            {
                for ( auto pIt = indicesIt->second.begin(); pIt != indicesIt->second.end(); ++pIt )
                {
                    f << *pIt << " ";
                } //...for each point index
                f << "\n\n";
            } //...if primitive has points
            else
                std::cerr << "[" << __func__ << "]: " << "primitive " << primitiveId << " does not have assigned points" << std::endl;
        } //...each primitive
        f << "# End of Primitives\n";

        f.close();
        std::cout << "[" << __func__ << "]: " << "wrote to " << outPath << std::endl;
        return EXIT_SUCCESS;
    } //...toGlobFit()
}

#endif // GF2_GLOBFIT_HPP
