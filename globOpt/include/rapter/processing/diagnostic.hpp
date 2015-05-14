#ifndef RAPTER_DIAGNOSTIC_HPP
#define RAPTER_DIAGNOSTIC_HPP

#include <string>
#include <map>
#include "Eigen/Sparse"
#include "rapter/simple_types.h"

namespace rapter
{

    template <typename _Scalar>
    class Diagnostic
    {
        public:
            typedef _Scalar                                      Scalar;
            typedef Eigen::SparseMatrix<_Scalar,Eigen::RowMajor> SparseMatrix;
            typedef std::map<LidT,std::string>                   NodeMapT;
            typedef std::pair<LidT,LidT>                         EdgeT;
            typedef std::vector< EdgeT >                         EdgeListT;
            typedef std::map< LidT, Eigen::Matrix<_Scalar,2,1> > NodePositionsT;

            Diagnostic( SparseMatrix const& qo, SparseMatrix const& Qo )
                : _qo( qo ), _Qo( Qo ) {}

            void setNodeName( LidT nodeId, std::string name )
            {
                _nodeNames[ nodeId ] = name;
            }

            template <class _Derived>
            void setNodePos( LidT nodeId, _Derived const& pos )
            {
                _nodePos[ nodeId ] = pos.template head(2).template cast<_Scalar>();
            }

            void addEdge( LidT v0, LidT v1 )
            {
                _edges.push_back( EdgeT(v0,v1) );
            }

            void draw( std::string path, bool show = false )
            {
                std::ofstream f;
                f.open( path.c_str() );
                if ( !f.is_open() ) { std::cerr << "could not open " << path << std::endl; return; }
                f << "graph {\n";

                char name[256];
                for ( typename NodeMapT::const_iterator it = _nodeNames.begin(); it != _nodeNames.end(); ++it )
                {
                    // first: nodeId
                    // second: nodeName
                    const LidT lid  = it ->first;

                    if ( _qo.coeff(lid,0) > _Scalar(0.) )
                    {
                        sprintf( name, "c%ld", lid );
                        f << "\t\"" << it->second << "\""
                          << " -- "
                          << "\"" << name << "\"\n";
                        f << "\t\"" << name << "\" [ label=\"" << _qo.coeff(lid,0) << "\"";
                        f << "];\n";
                        if ( _nodePos.find(lid) != _nodePos.end() )
                        {
                            f << "\t\"" << it->second << "\"[ pos=\"" << _nodePos[lid](0) * 10. << ","
                                                                      << _nodePos[lid](1) * 10. << "!\"";
                            f << "];\n";
                        }
                    }
                    for ( typename NodeMapT::const_iterator it2 = _nodeNames.begin(); it2 != _nodeNames.end(); ++it2 )
                    {
                        // first: nodeId
                        // second: nodeName
                        const LidT lid1 = it2->first;
                        if ( lid1 <= lid ) continue;

                        if ( (_Qo.coeff(lid,lid1) > _Scalar(0.)) || (_Qo.coeff(lid1,lid) > _Scalar(0.)) )
                        {
                            f << "\t\"" << it->second << "\""
                              << " -- "
                              << "\"" << it2->second << "\""
                              << " [label=\"" << _Qo.coeff(lid1,lid) + _Qo.coeff(lid,lid1) << "\"]\n";
                        }
                    }
                }

                for ( EdgeListT::const_iterator eIt = _edges.begin(); eIt != _edges.end(); ++eIt )
                {
                    if ( (_nodeNames.find(eIt->first) != _nodeNames.end()) &&
                         (_nodeNames.find(eIt->second) != _nodeNames.end()) )
                    f << "\t\"" << _nodeNames.at(eIt->first) << "\""
                      << " -- "
                      << "\"" << _nodeNames.at(eIt->second) << "\"\n";
                }

                f << "}\n";
                f.close();

                if ( show )
                {
                    char cmd[2048];
                    sprintf( cmd, "(dot -Tpng -o%s.png %s &)", path.c_str(), path.c_str() );
                    int sysErr = system( cmd );
                    if ( sysErr ) std::cout << "sysErr: " << sysErr << std::endl;
                }
            }

        protected:
            SparseMatrix    _qo, _Qo;
            NodeMapT        _nodeNames;
            EdgeListT       _edges;
            NodePositionsT  _nodePos;

//            inline void _addEdge( std::ofstream &f )
//            {
//                f
//            }
    };

} //...ns rapter

#endif // RAPTER_DIAGNOSTIC_HPP
