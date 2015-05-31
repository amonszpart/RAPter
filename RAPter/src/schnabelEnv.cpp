#ifndef __RAPTER_SCHNABELENV_HPP__
#define __RAPTER_SCHNABELENV_HPP__

#include "schnabelEnv.h"

#include "pcl/visualization/pcl_visualizer.h"
//#include "pcltools/util.hpp" // smartgeometry
//#include "AMUtil2.h"

#include "rapter/primitives/planePrimitive.h"

// --- schnabel07 ---
#include "PointCloud.h"
#include "RansacShapeDetector.h"
#include "MiscLib/RefCountPtr.h"
#include "PlanePrimitiveShapeConstructor.h"
#include "PlanePrimitiveShape.h"
// --- END schnabel07 ---

namespace rapter
{
    template <class PclCloudT, typename PrimitiveT, /*class PidGidT,*/ class PointContainerT >
    int SchnabelEnv::run( std::vector<PrimitiveT>          &planes
                          , PointContainerT &outPoints //PidGidT &pidGid
                          , PointContainerT const& points
                          , typename PclCloudT::Ptr           &inCloud
                          , float scale
                          , int                                     min_support_arg
                          , int show
                          , bool extrude2D
                          , int pointMultiplier
                          )
    {
        typedef typename PointContainerT::value_type PointPrimitiveT;
        typedef typename PclCloudT::PointType PT;

        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
        {
            vptr->setBackgroundColor( .4, .8, .4 );
        }

        // 2D extrude
        typename PclCloudT::Ptr cloud;
        {
            if ( extrude2D )
            {
                cloud.reset( new PclCloudT() );
                for (PidT pid = 0; pid < inCloud->size(); ++pid )
                {
                    for ( int i = 0; i != pointMultiplier; ++i )
                    {
                        cloud->push_back( inCloud->at( pid ) );
                        cloud->back().z += rand() / float(RAND_MAX);
                    }
                }
            }
            else
                cloud = inCloud;
        }
        std::cout << "[" << __func__ << "]: " << "finished cloud copy" << std::endl;

        vptr->addPointCloud<PT>( cloud, "xtruded" );

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
        //schnabel::Point pnts[ cloud->size() ];
        schnabel::Point* pnts = new schnabel::Point[ cloud->size() ];
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
        MiscLib::Vector< int > outShapeIndex; // points->planes assignments output
        std::vector< MiscLib::Vector< size_t > > outIndices;
        {
            std::cout << "starting..." << std::endl;
            int ret = rsd.Detect( pc, 0, pc.size(), &shapes, &outShapeIndex, &outIndices );
            std::cout << "detect returned " << ret << std::endl;
            std::cout << "shapes.size: " << shapes.size() << std::endl;
        }

        outPoints.clear(); outPoints.reserve( outShapeIndex.size() );
        size_t eidx = pc.size();
        for ( LidT shapeId = 0; shapeId != outIndices.size(); ++shapeId )
        {
            size_t bidx = eidx - shapes[shapeId].second;
            for ( PidT pidId = bidx; pidId < eidx; ++pidId )
            {
                PointPrimitiveT pnt( Eigen::Map<const Eigen::Vector3f>( pc.at(pidId).pos.getValue(), 3 )
                                   , Eigen::Map<const Eigen::Vector3f>( pc.at(pidId).normal        , 3 ) );
                pnt.setTag( PointPrimitiveT::TAGS::PID, outPoints.size() );
                pnt.setTag( PointPrimitiveT::TAGS::GID, shapeId );
                //schnabel::Point( Vec3f(cloud->at(i).x, cloud->at(i).y, cloud->at(i).z) );
                outPoints.push_back( pnt );
            }
            eidx = bidx;
        }
#if 0
            std::cout << "shape->size: " << shapes[shapeId]->second << ", outIndices[" << shapeId<< "]: " << outIndices[shapeId].size() << std::endl;
            for ( PidT pidId = 0; pidId != outIndices[shapeId].size(); ++pidId, ++curr )
            {
                PidT pid = outIndices[shapeId][pidId];
                ///std::cout << "outsh[" << pid << "]: " << outShapeIndex[ pid ] << std::endl;
                //pidGid[ pid ] = /* planeId, converted to GID later: */ outShapeIndex[ pid ];
                //            std::cout << "pc.at" << pid << ": " << pc.at( pid ).pos.getValue()[0]
                //                                        << ", " << pc.at(pid).normal << std::endl;
                PointPrimitiveT pnt( Eigen::Map<const Eigen::Vector3f>( pc.at(curr).pos.getValue(), 3 )
                                   , Eigen::Map<const Eigen::Vector3f>( pc.at(curr).normal        , 3 ) );
                //PointPrimitiveT pnt = points.at( pid );
                pnt.setTag( PointPrimitiveT::TAGS::PID, outPoints.size() );
                pnt.setTag( PointPrimitiveT::TAGS::GID, outIndices.size() - shapeId - 1 );
                //schnabel::Point( Vec3f(cloud->at(i).x, cloud->at(i).y, cloud->at(i).z) );
                outPoints.push_back( pnt );
            }
        }
#endif

        // process output
        for ( size_t i = 0; i != shapes.size(); ++i )
        {
            // sort primitive type
            switch ( shapes[i].first->Identifier() )
            {
                case schnabel::PrimitiveShape::PLANE_PRIMITIVE_SHAPE:
                {
                    MiscLib::RefCountPtr<schnabel::PlanePrimitiveShape> planePS = static_cast<schnabel::PlanePrimitiveShape*>( shapes[i].first->Clone() );
                    std::cout << "shapes[" << i << "].second: " << shapes[i].second << std::endl;
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
                    //planes.back()()(2) *= 2.;
                    //std::cout << "dist to orig: " << planes.back().getDistance( Eigen::Vector3f::Zero() ) << std::endl;

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
#if 0
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
            points[pid].setTag( PointPrimitiveT::TAGS::GID, min_gid );
        }
        std::cout << "finishing assignment" << std::endl;
#endif

        if ( pnts ) delete [] pnts;
        return EXIT_SUCCESS;
    }

} // ns am


#endif // __RAPTER_SCHNABELENV_HPP__
