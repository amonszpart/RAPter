//=======================================================================
// Copyright 2001 Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee,
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

#ifndef GF2_MST_HPP
#define GF2_MST_HPP

#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

namespace rapter
{

namespace mst
{
    typedef std::pair< int, int > Edge;
    template <typename _Scalar>
    class MST
    {
        public:
            //addEdge()
    };

template <typename _Scalar>
int mstMain()
{
    using namespace boost;

    typedef adjacency_list < vecS
                           , vecS
                           , undirectedS
                           , property<vertex_distance_t, _Scalar>
                           , property < edge_weight_t, _Scalar >
                           > Graph;


    const int num_nodes = 5;
    Edge edges[] = { Edge(0, 2), Edge(1, 3), Edge(1, 4), Edge(2, 1), Edge(2, 3),
                  Edge(3, 4), Edge(4, 0)
                };
    _Scalar weights[] = { 1, 1, 2, 7, 3, 1, 1 };
    int num_edges = sizeof(edges) / sizeof(Edge);
    Graph g(edges, edges + num_edges, weights, num_nodes);
    typename property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
    std::vector < typename graph_traits < Graph >::vertex_descriptor > p( num_vertices(g) );

    prim_minimum_spanning_tree(g, &p[0]);
    std::fstream dotfile("mst.dot");
    if ( !dotfile.is_open() )
    {
        std::cerr << "[" << __func__ << "]: " << "could not open mst.dot" << std::endl;
        return EXIT_FAILURE;
    }

    dotfile << "digraph {\nedge [dir=none]\n";
    // Original
    for ( std::size_t i = 0; i != num_edges; ++i )
    {
        dotfile << "V" << edges[i].first << " -> " << "V" << edges[i].second << std::endl;
    }
    // MST
    for (std::size_t i = 0; i != p.size(); ++i)
        if (p[i] != i)
        {
            std::cout << "parent[" << i << "] = " << p[i] << std::endl;
            dotfile << "v" << p[i] << " -> " << "v" << i << std::endl;
        }
        else
            std::cout << "parent[" << i << "] = no parent" << std::endl;

    dotfile << "}\n";
    dotfile.close();
    std::system("dot -Tpng mst.dot -o mst.png && eog mst.png");

    return EXIT_SUCCESS;
}

} //...ns mst

} //...ns GF2

#endif // GF2_MST_HPP
