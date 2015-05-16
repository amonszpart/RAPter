#ifndef GO_TRIANGLESFROMOBJ_HPP
#define GO_TRIANGLESFROMOBJ_HPP

#include "pcl/io/obj_io.h"
#include "pcl/PolygonMesh.h"

namespace rapter
{

    template <class _PclCloud, typename _TriangleT>
    inline int getTriangleFromFace( _TriangleT & triangle, _PclCloud const& cloud, ::pcl::Vertices const& face )
    {
        typedef typename _TriangleT::Scalar Scalar;
        typedef typename _TriangleT::VectorT Vector;

        triangle = _TriangleT( cloud.at(face.vertices[0]).template getVector3fMap().template cast<Scalar>(),
                               cloud.at(face.vertices[1]).template getVector3fMap().template cast<Scalar>(),
                               cloud.at(face.vertices[2]).template getVector3fMap().template cast<Scalar>() );
        return EXIT_SUCCESS;
    } //...getFace

    template <class _TrianglesContainer, typename Scalar>
    inline int getTrianglesFromObj( _TrianglesContainer & triangles, std::string const& meshPath, Scalar minPlaneEdge )
    {
        typedef typename _TrianglesContainer::value_type Triangle;
        bool filterBySize = minPlaneEdge > 0.;

        // get mesh
        pcl::PolygonMesh mesh;
        pcl::io::loadOBJFile( meshPath, mesh );

        // parse input mesh
        size_t cnt = 0;
        if ( meshPath.find("obj") < 0 )
        {
            std::cout << "--assign " << meshPath << " has to be obj and contain triangles!" << std::endl;
            return EXIT_FAILURE;
        }
        else
        {
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
                getTriangleFromFace( tri, polyCloud, *it );
                std::vector<Scalar> sideLengths = tri.getSideLengths();
                if ( filterBySize && *std::min_element(sideLengths.begin(), sideLengths.end()) < minPlaneEdge )
                {
                    //std::cout << "filtering triangle, since minedge: " << *std::min_element(sideLengths.begin(), sideLengths.end()) << " < " << minPlaneEdge << std::endl;
                    ++cnt;
                }
                else
                    triangles.push_back( tri );

                //triangles.push_back( Triangle(it->vertices[0], it->verticess[1], it->vertices[2]) );
            } //...for polygons
        } //...read obj

        std::cout << "have " << triangles.size() << " triangles" << ", filtered " << cnt << " = " << float(cnt)/mesh.polygons.size() << std::endl;
        return EXIT_SUCCESS;
    }

} //...ns globopt

#endif // GO_TRIANGLESFROMOBJ_HPP
