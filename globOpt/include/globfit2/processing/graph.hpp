#ifndef GF2_GRAPH_HPP
#define GF2_GRAPH_HPP

#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace GF2
{
    /*! \brief Namespace graph to make usage of boost graph class wrapper \ref Graph easier.
     */
    namespace graph
    {
        /*! \brief tuple<3> to temporarily store graph edges until we know how large the graph will be
         *  \tparam _Scalar EdgeWeight type.
         */
        template <typename _Scalar>
        struct EdgeT
        {
            int _v0, _v1;
            _Scalar _w;

            EdgeT( int v0, int v1, _Scalar w ) : _v0(v0), _v1(v1), _w(w) {}

            bool operator<( EdgeT<_Scalar> const& other ) const
            {
                if ( _v0 < other._v0 ) return true;
                else if ( _v1 < other._v1 ) return true;
                else return _w < other._w;
            }
        };

        /*! \brief Class instead of "typedef std::set< graph::EdgeT<_Scalar> >".
         *         Records maximum vertex Id in the edgelist.
         * \tparam _Scalar EdgeWeight type.
         */
        template <typename _Scalar>
        class EdgeListT : public std::set<graph::EdgeT<_Scalar> >
        {
            typedef std::set<graph::EdgeT<_Scalar> > ParentT;

            public:
                /*! \brief Overload constructor to initialize vertex id field.
                 */
                EdgeListT() : ParentT(), _maxVertexId( UidT(0) ) {}

                /*! \brief Overload set insert to keep track of maximum node id in the edges list.
                 */
                std::pair<typename ParentT::iterator, bool>
                insert( typename ParentT::value_type&& __x )
                {
                    if ( __x._v0 > _maxVertexId ) _maxVertexId = __x._v0;
                    if ( __x._v1 > _maxVertexId ) _maxVertexId = __x._v1;

                    return ParentT::insert( __x );
                }

                /*! \brief Get maximum vertex id added up till now.
                 *  \warning Does not clear upon clear!! \todo Clear vertex id upon clear;
                 */
                inline UidT getMaxVertexId() const { return _maxVertexId; }

            protected:
                UidT _maxVertexId;
        };

    } //...ns graph

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

            typedef std::vector<int> ComponentListT;
            typedef std::map<int,size_t> ComponentSizesT;

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

            Graph( graph::EdgeListT<_Scalar> const& edgesList )
                : Graph( edgesList.getMaxVertexId() )
            {
                for ( auto it = edgesList.begin(); it != edgesList.end(); ++it )
                    this->addEdge( it->_v0, it->_v1, /* not used right now: */ it->_w );
            }

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
            inline int getComponents( std::vector<int> &components, ComponentSizesT *counts = NULL )
            {
                if ( counts ) counts->clear();

                components.resize( boost::num_vertices(_undigraph) );
                int num = boost::connected_components(_undigraph, &components[0]);

                std::vector<int>::size_type i;
                for (i = 0; i != components.size(); ++i)
                {
                    if ( counts ) (*counts)[ components[i] ]++;
                }
                std::cout << std::endl;

                return num;
            }

            inline int draw( std::string const& out_path, bool show = false )
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

                            //std::cout << ss.str();
                            f << ss.str();
                        }
                    }
                }

                f << "}\n";
                f.close();

                if ( show )
                {
                    char cmd[2048];
                    sprintf( cmd, "dot -Tpng -o %s.png %s && (eog %s.png &)", out_path.c_str(), out_path.c_str(), out_path.c_str() );
                    std::cout << "[" << __func__ << "]: " << cmd << std::endl;
                    system( cmd );
                }
            }

            typedef std::map< int, std::vector<UidT> > ClustersT; // [ cluster0: [v0, v10,...], cluster1: [v3, v5, ...], ... ]

            inline void getClusters( ClustersT &clusters, int sizeLimit = 2 )
            {
                ComponentListT      components;
                ComponentSizesT     compSizes;
                this->getComponents( components, &compSizes );

                for ( UidT uId = 0; uId != components.size(); ++uId )
                {
                    if ( compSizes[ components[uId] ] < sizeLimit )
                        continue;

                    // store
                    clusters[ components[uId] ].push_back( uId );
                }
            } //...getClusters()

            inline void showClusters( ClustersT const& clusters,
                                      std::string const& path, bool show = false )
            {
                std::ofstream f;
                f.open( path.c_str() );
                f << "graph {\n";

                for ( ClustersT::const_iterator it = clusters.begin(); it != clusters.end(); ++it )
                {
                    // it->first: clusterId
                    // it->second: vector<UidT>
                    for ( typename ClustersT::mapped_type::const_iterator nodeIt = it->second.begin();
                          nodeIt != it->second.end(); ++nodeIt )
                    {
                        f << *nodeIt << " -- c" << /* clusterId */ it->first << std::endl;
                    }
                }

                f << "}" << std::endl;
                f.close();

                if ( show )
                {
                    char cmd[2048];
                    sprintf( cmd, "dot -Tpng -o %s.png %s && (eog %s.png &)", path.c_str(), path.c_str(), path.c_str() );
                    std::cout << "[" << __func__ << "]: " << cmd << std::endl;
                    system( cmd );
                }
            }

        protected:
            _UndirectedGraph _undigraph;
            std::map<int,std::string> _names;
    };

} //...ns GF2

#endif // GF2_GRAPH_HPP
