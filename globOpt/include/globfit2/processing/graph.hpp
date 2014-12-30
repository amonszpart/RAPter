#ifndef GF2_GRAPH_HPP
#define GF2_GRAPH_HPP

#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace GF2
{
    template <typename _Scalar>
    struct MyGraphConfig
    {
        typedef boost::property<boost::edge_weight_t, _Scalar> WeightT;
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, WeightT> UndirectedGraph;
    };

    template <typename _Scalar, typename _UndirectedGraph /* = MyGraphConfig::UndirectedGraph*/>
    class Graph
    {
        public:
            typedef typename boost::graph_traits< _UndirectedGraph >::vertex_descriptor VertexT;
            typedef typename boost::graph_traits< _UndirectedGraph >::edge_descriptor EdgeT;
            typedef typename boost::graph_traits< _UndirectedGraph >::out_edge_iterator OutEdgeIteratorT;

            static inline int testGraph()
            {
                const int V = 3;
                  _UndirectedGraph undigraph(V);
                  VertexT zero, one, two;
                  OutEdgeIteratorT out, out_end;
                  typename boost::graph_traits< _UndirectedGraph >::in_edge_iterator in, in_end;

                  zero = boost::vertex(0, undigraph);
                  one = boost::vertex(1, undigraph);
                  two = boost::vertex(2, undigraph);
                  boost::add_edge(zero, one, undigraph);
                  boost::add_edge(zero, two, undigraph);
                  boost::add_edge(one, two, undigraph);

                  std::cout << "out_edges(0):";
                  for (boost::tie(out, out_end) = out_edges(zero, undigraph); out != out_end; ++out)
                    std::cout << ' ' << *out;
                  std::cout << std::endl << "in_edges(0):";
                  for (boost::tie(in, in_end) = in_edges(zero, undigraph); in != in_end; ++in)
                    std::cout << ' ' << *in;
                  std::cout << std::endl;
            }

            Graph( const int vertex_count )
                : _undigraph( vertex_count )
            {}

            inline std::pair<EdgeT, bool>
            addEdge( const int v0, const int v1, _Scalar const weight )
            {
                //throw new std::runtime_exception("todo" );
                return boost::add_edge( v0, v1, _undigraph );
                //return 0;
            }

            inline std::string
            addVertexName( const int v, std::string const& name )
            {
                if ( _names.find(v) != _names.end() )
                    std::cerr << "[" << __func__ << "]: " << "adding name to already named node!" << std::endl;

                _names[v] = name;

                return _names[v];
            } //...addVertexName

            /*! \brief Gets the graph's connected components.
             *  \param[in] counts Entries: <component_id, component_size>
             */
            inline int getComponents( std::vector<int> &components, std::map<int, int> *counts = NULL )
            {
                if ( counts ) counts->clear();

                components.resize( boost::num_vertices(_undigraph) );
                int num = boost::connected_components(_undigraph, &components[0]);

                std::vector<int>::size_type i;
                std::cout << "Total number of components: " << num << std::endl;
                for (i = 0; i != components.size(); ++i)
                {
                    std::cout << "Vertex " << i <<" is in component " << components[i] << std::endl;
                    if ( counts ) (*counts)[ components[i] ]++;
                }
                std::cout << std::endl;

                return num;
            }

            inline int draw( std::string const& out_path )
            {
                std::ofstream f;
                f.open( out_path );
                if ( !f.is_open() ) { std::cerr <<  "[" << __func__ << "]: " << "could not open " << out_path  << std::endl; return 1; }
                f << "graph {\n";

                OutEdgeIteratorT out,out_end;
                for ( int vid = 0; vid != boost::num_vertices(_undigraph); ++vid )
                {
                    for ( boost::tie(out, out_end) = boost::out_edges(boost::vertex(vid,_undigraph), _undigraph); out != out_end; ++out )
                    {
                        int v0 = boost::source(*out,_undigraph), v1 = boost::target(*out,_undigraph);

                        if ( v0 < v1 )
                        {
                            auto nameIt0 = _names.find( v0 );
                            auto nameIt1 = _names.find( v1 );
                            std::stringstream ss;
                            if ( nameIt0 != _names.end() )
                                 ss << nameIt0->second;
                            else
                                 ss << v0;
                            ss << " -- ";
                            if ( nameIt1 != _names.end() )
                                 ss << nameIt1->second;
                            else
                                 ss << v1;
                            ss << "\n";

                            std::cout << ss.str();
                            f << ss.str();
                        }
                    }
                }

                f << "}\n";
                f.close();
            }

        protected:
            _UndirectedGraph _undigraph;
            std::map<int,std::string> _names;
    };

} //...ns GF2

#endif // GF2_GRAPH_HPP
