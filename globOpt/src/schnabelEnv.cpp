#include "schnabelEnv.h"

#include "pcl/visualization/pcl_visualizer.h"
#include "pcltools/util.hpp" // smartgeometry
//#include "AMUtil2.h"

#include "globfit2/primitives/planePrimitive.h"

// --- schnabel07 ---
#include "PointCloud.h"
#include "RansacShapeDetector.h"
#include "MiscLib/RefCountPtr.h"
#include "PlanePrimitiveShapeConstructor.h"
#include "PlanePrimitiveShape.h"
// --- END schnabel07 ---

namespace GF2
{
    template <typename PrimitiveT>
    int SchnabelEnv::run( std::vector<PrimitiveT>          &planes
                        , pcl::PointCloud<GF2::MyPoint>::Ptr      cloud
                        , float scale
                        , int                                     min_support_arg
                          , int show )
    {
        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
        {
            vptr->setBackgroundColor( .4, .8, .4 );
        }

        // read
        pcl::PointCloud<pcl::Normal>::Ptr normals;
        smartgeometry::calculateNormals<MyPoint>( normals
                                                  , cloud
                                                  , /*         indices: */ nullptr
                                                  , /*       normalize: */ true
                                                  , /* K_neighbourhood: */ 20 );

        // show
        vptr->addPointCloud<GF2::MyPoint>( cloud );
        vptr->addCoordinateSystem( 0.1, "coordsys", 0 );
        vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4.f );

        // Points
        Point pnts[ cloud->size() ];
        for ( size_t i = 0; i != cloud->size(); ++i )
        {
            pnts[i] = Point(   Vec3f( cloud  ->at(i).       x, cloud  ->at(i).       y, cloud  ->at(i).       z) );
        }
        PointCloud pc( pnts, cloud->size() );

        // Normals
        pc.calcNormals( scale );
        for ( size_t i = 0; i != cloud->size(); ++i )
        {
            // error check
            if ( (cloud->at(i).x != pc.at(i)[0]) || (cloud->at(i).y != pc.at(i)[1]) || (cloud->at(i).z != pc.at(i)[2]) )
            {
                std::cerr << cloud->at(i).x << " != " << pc.at(i)[0] << "\t || \t"
                          << cloud->at(i).y << " != " << pc.at(i)[1] << "\t || \t"
                          << cloud->at(i).z << " != " << pc.at(i)[2] << std::endl;
            }

            // debug pcl normals
            //                std::cout << "normal" << i << ": "
            //                          << normals->at(i).normal[0] << ", "
            //                          << normals->at(i).normal[1] << ", "
            //                          << normals->at(i).normal[2] << "\t =? \t "
            //                          << pc.at(i).normal.print() << std::endl;

            // copy for display: schnabel -> pcl
            normals->at(i).normal[0] = pc.at(i).normal[0];
            normals->at(i).normal[1] = pc.at(i).normal[1];
            normals->at(i).normal[2] = pc.at(i).normal[2];
            Eigen::Map<Eigen::Vector3f>( normals->at(i).normal ).normalize();

            //                std::cout << "normal" << i << ": "
            //                          << normals->at(i).normal[0] << ", "
            //                          << normals->at(i).normal[1] << ", "
            //                          << normals->at(i).normal[2] << "\t =? \t "
            //                          << pc.at(i).normal.print() << std::endl;

        }

        // show normals
        vptr->addPointCloudNormals<GF2::MyPoint,pcl::Normal>( cloud, normals, 10, 0.15, "cloud_normals2", 0 );
        vptr->spin();

        // options (schnabel)
        RansacShapeDetector::Options opt;
        opt.m_minSupport = min_support_arg;
        RansacShapeDetector rsd( opt );
        rsd.Add( new PlanePrimitiveShapeConstructor() );

        // work
        MiscLib::Vector< std::pair< MiscLib::RefCountPtr< PrimitiveShape >, size_t > > shapes; // output
        {
            std::cout << "starting..." << std::endl;
            int ret = rsd.Detect( pc, 0, pc.size(), &shapes );
            std::cout << "detect returned " << ret << std::endl;
            std::cout << "shapes.size: " << shapes.size() << std::endl;
        }

        // process output
        for ( size_t i = 0; i != shapes.size(); ++i )
        {
            // sort primitive type
            switch ( shapes[i].first->Identifier() )
            {
                case PrimitiveShape::PLANE_PRIMITIVE_SHAPE:
                {
                    MiscLib::RefCountPtr<PlanePrimitiveShape> planePS = static_cast<PlanePrimitiveShape*>( shapes[i].first->Clone() );
                    Vec3f p,n;
                    planePS->Normal( p, &n );
                    float dist = -1.f * planePS->Internal().SignedDistToOrigin();
                    //std::cout << "\tnormal: " << n.print() << std::endl;
                    //std::cout << "\tdist: " << dist << std::endl;
                    //std::cout << "\tpos: " << planePS->Internal().getPosition().print() << std::endl;

                    pcl::PointXYZ centroid( planePS->Internal().getPosition()[0],
                                            planePS->Internal().getPosition()[1],
                                            planePS->Internal().getPosition()[2] );

                    //planes.emplace_back( PlanePrimitive((Eigen::Vector4f() << n[0],n[1],n[2],dist).finished()) );
                    Eigen::Vector3f normal; normal << n[0],n[1],n[2];
                    //planes.emplace_back( PlanePrimitive(Eigen::Vector3f::Zero() + normal * dist, normal) );
                    planes.emplace_back( PrimitiveT(centroid.getVector3fMap(), normal) );

                    char name[255];
                    sprintf(name, "plane%d",i);
                    vptr->addPlane( *(planes.back().modelCoefficients()),
                                    centroid.x,centroid.y,centroid.z,
                                    name, 0);

                    //                            vptr->addArrow( centroid,
                    //                                            pcl::PointXYZ(0.f,0.f,0.f),
                    //                                            false,
                    //                                            1.f,0.f,0.f,
                    //                                            am::util::sprintf("planepos%d",i),
                    //                                            0 );
                    //std::cout << "added plane: " << planes.back().coeffs().transpose() << std::endl;
                    //vptr->spin();
                    break;
                }

                default:
                {
                    std::cerr << __func__ << "unknown shape identifier..." << std::endl;
                    break;
                }
            }
        }

        vptr->spin();

        return EXIT_SUCCESS;
    }

} // ns am
