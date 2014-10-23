#ifndef GF2_SCHNABELENV_HPP_
#define GF2_SCHNABELENV_HPP_

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
    template <class PclCloudT, typename PrimitiveT, class PidGidT, class PointContainerT >
    int SchnabelEnv::run( std::vector<PrimitiveT>          &planes
                          , PidGidT &pidGid
                          , PointContainerT &points
                          , typename PclCloudT::Ptr           &cloud
                          , float scale
                          , int                                     min_support_arg
                          , int show )
    {
        typedef typename PointContainerT::value_type PointPrimitiveT;
        typedef typename PclCloudT::PointType PT;

        for ( size_t i = 0; i != std::min(10UL,points.size()); ++i )
        {
            std::cout << "\tpoints[" << i << "]: " << points[i].toString() << "\n\tcloud[" << i << "]: " << cloud->at(i).getVector3fMap().transpose() << std::endl;
        }

        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
        {
            vptr->setBackgroundColor( .4, .8, .4 );
        }

        // read
        pcl::PointCloud<pcl::Normal>::Ptr normals( new pcl::PointCloud<pcl::Normal> );
        for ( int i = 0; i != cloud->size(); ++i )
        {
            pcl::Normal nrm;
            nrm.normal_x = cloud->at(i).normal_x;
            nrm.normal_y = cloud->at(i).normal_y;
            nrm.normal_z = cloud->at(i).normal_z;
            normals->push_back( nrm );
        }
//        smartgeometry::calculateNormals<MyPoint>( normals
//                                                  , cloud
//                                                  , /*         indices: */ nullptr
//                                                  , /*       normalize: */ true
//                                                  , /* K_neighbourhood: */ 20 );

        // show
        vptr->addPointCloud<PT>( cloud );
        vptr->addCoordinateSystem( 0.1, "coordsys", 0 );
        vptr->setPointCloudRenderingProperties( pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4.f );

        // Points
        schnabel::Point pnts[ cloud->size() ];
        for ( size_t i = 0; i != cloud->size(); ++i )
        {
            pnts[i] = schnabel::Point( Vec3f(cloud->at(i).x, cloud->at(i).y, cloud->at(i).z) );
        }
        schnabel::PointCloud pc( pnts, cloud->size() );

        // Normals
        pc.calcNormals( scale );
        for ( size_t i = 0; i != cloud->size(); ++i )
        {
            // error check
//            if ( (cloud->at(i).x != pc.at(i)[0]) || (cloud->at(i).y != pc.at(i)[1]) || (cloud->at(i).z != pc.at(i)[2]) )
//            {
//                std::cerr << cloud->at(i).x << " != " << pc.at(i)[0] << "\t || \t"
//                          << cloud->at(i).y << " != " << pc.at(i)[1] << "\t || \t"
//                          << cloud->at(i).z << " != " << pc.at(i)[2] << std::endl;
//            }

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
        if ( show )
        {
            //vptr->addPointCloudNormals<PT,pcl::Normal>( cloud, normals, 10, 0.15, "cloud_normals2", 0 );

            //vptr->spin();
        }

        // options (schnabel)
        schnabel::RansacShapeDetector::Options opt;
        opt.m_minSupport = min_support_arg;
        schnabel::RansacShapeDetector rsd( opt );
        rsd.Add( new schnabel::PlanePrimitiveShapeConstructor() );

        // work
        MiscLib::Vector< std::pair< MiscLib::RefCountPtr< schnabel::PrimitiveShape >, size_t > > shapes; // output
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
                case schnabel::PrimitiveShape::PLANE_PRIMITIVE_SHAPE:
                {
                    MiscLib::RefCountPtr<schnabel::PlanePrimitiveShape> planePS = static_cast<schnabel::PlanePrimitiveShape*>( shapes[i].first->Clone() );
                    Vec3f p,n;
                    planePS->Normal( p, &n );
                    float dist = -1.f * planePS->Internal().SignedDistToOrigin();
                    std::cout << "\tnormal: " << n.print() << std::endl;
                    std::cout << "\tdist: " << dist << std::endl;
                    std::cout << "\tpos: " << planePS->Internal().getPosition().print() << std::endl;


                    pcl::PointXYZ centroid; centroid.getVector3fMap() <<  planePS->Internal().getPosition()[0],
                                            planePS->Internal().getPosition()[1],
                                            planePS->Internal().getPosition()[2];

                    //planes.emplace_back( PlanePrimitive((Eigen::Vector4f() << n[0],n[1],n[2],dist).finished()) );
                    Eigen::Vector3f normal; normal << n[0],n[1],n[2];
                    //planes.emplace_back( PlanePrimitive(Eigen::Vector3f::Zero() + normal * dist, normal) );
                    planes.emplace_back( PrimitiveT(centroid.getVector3fMap(), normal) );
                    planes.back()()(2) *= 2.;
                    std::cout << "dist to orig: " << planes.back().getDistance( Eigen::Vector3f::Zero() ) << std::endl;

                    if ( show && i < 6)
                    {
                        char name[255];
                        sprintf(name, "plane%lu",i);
                        vptr->addPlane( *(planes.back().modelCoefficients()),
                                        centroid.x,centroid.y,centroid.z,
                                        name, 0);
                        vptr->addArrow( centroid, pcl::PointXYZ(0.5,0.5,0.5), 0., 1., .5, false, name + std::string("arrow"), 0 );
                    }

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

        if ( show )
            vptr->spin();
        else
            vptr->spinOnce();

        // assign points
#if 1
        std::cout << "starting assignment" << std::endl; fflush( stdout );
        for ( int pid = 0; pid != points.size(); ++pid )
        {
            float min_dist = FLT_MAX, tmp; int min_gid = 0;
            for ( int gid = 0; gid != planes.size(); ++gid )
            {
                tmp = std::abs( planes[gid].getDistance(points[pid].template pos()) );
                if ( tmp < 0.f )
                    std::cerr << "asdf: " << tmp << std::endl;
                if ( (tmp < min_dist) && (tmp < scale) )
                {
                    min_dist = tmp;
                    min_gid = gid;
                }
            }
            pidGid[ pid ] = min_gid;
            points[pid].setTag( PointPrimitiveT::GID, min_gid );
        }
        std::cout << "finishing assignment" << std::endl;
#endif

        return EXIT_SUCCESS;
    }

} // ns am


#endif //GF2_SCHNABELENV_HPP_
