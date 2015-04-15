#ifndef GO_TRIANGLESFROMOBJ_HPP
#define GO_TRIANGLESFROMOBJ_HPP

#include "pcl/io/obj_io.h"
#include "pcl/PolygonMesh.h"

namespace globopt
{

    template <class _PclCloud, typename _TriangleT>
    inline int getFace( _TriangleT & triangle, _PclCloud const& cloud, ::pcl::Vertices const& face )
    {
        typedef typename _TriangleT::Scalar Scalar;
        typedef typename _TriangleT::VectorT Vector;

        triangle = _TriangleT( cloud.at(face.vertices[0]).template getVector3fMap().template cast<Scalar>(),
                cloud.at(face.vertices[1]).template getVector3fMap().template cast<Scalar>(),
                cloud.at(face.vertices[2]).template getVector3fMap().template cast<Scalar>() );
        return EXIT_SUCCESS;
    } //...getFace

    template <class _TrianglesContainer>
    inline int getTrianglesFromObj( _TrianglesContainer & triangles, std::string const& meshPath )
    {
        typedef typename _TrianglesContainer::value_type Triangle;

        // parse input mesh
        if ( meshPath.find("obj") < 0 )
        {
            std::cout << "--assign " << meshPath << " has to be obj and contain triangles!" << std::endl;
            return EXIT_FAILURE;
        }
        else
        {
            // get mesh
            pcl::PolygonMesh mesh;
            pcl::io::loadOBJFile( meshPath, mesh );
            // get mesh poins
            pcl::PointCloud<pcl::PointXYZRGB> polyCloud;
            pcl::fromPCLPointCloud2( mesh.cloud, polyCloud );

            // create Triangles from mesh faces
            triangles.reserve( mesh.polygons.size() );
            for ( auto it = mesh.polygons.begin(); it != mesh.polygons.end(); ++it )
            {
                if ( it->vertices.size() != 3 )
                {
                    std::cerr << "[" << __func__ << "]: " << "found a face with " << it->vertices.size() << " vertices, can't handle it!" << std::endl;
                    continue;
                }

                Triangle tri;
                getFace( tri, polyCloud, *it );
                triangles.push_back( tri );
                //triangles.push_back( Triangle(it->vertices[0], it->verticess[1], it->vertices[2]) );
            } //...for polygons
        } //...read obj

        std::cout << "have " << triangles.size() << " triangles" << std::endl;
        return EXIT_SUCCESS;
    }

} //...ns globopt

#endif // GO_TRIANGLESFROMOBJ_HPP
