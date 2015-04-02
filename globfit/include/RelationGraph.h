#ifndef RelationGraph_H
#define RelationGraph_H

#include <set>
#include <vector>
#include "boost/graph/adjacency_list.hpp"

#include "RelationEdge.h"

using namespace boost;

struct EdgeProp {
    RelationEdge::RelationEdgeType edgeType;
    double score;
    EdgeProp(RelationEdge::RelationEdgeType initEdgeType=RelationEdge::RET_PARALLEL, double initScore=1.0):edgeType(initEdgeType), score(initScore){}
};

typedef adjacency_list<setS, vecS, undirectedS, RelationVertex, EdgeProp> Graph;
typedef graph_traits<Graph>::vertex_descriptor   GraphVertex;
typedef graph_traits<Graph>::vertices_size_type  GraphVertexIdx;
typedef graph_traits<Graph>::edge_descriptor     GraphEdge;
typedef graph_traits<Graph>::vertex_iterator     GraphVertexItr;
typedef graph_traits<Graph>::edge_iterator       GraphEdgeItr;
typedef graph_traits<Graph>::out_edge_iterator   GraphOutEdgeItr;
typedef graph_traits<Graph>::adjacency_iterator  GraphAdjacencyItr;
typedef graph_traits<Graph>::out_edge_iterator   GraphOutEdgeItr;

void biconnectDecompose(std::vector<RelationEdge>& vecRelationEdge, Graph& g);

void reduceParaOrthEdges(const std::vector<Primitive*>& vecPrimitive, double angleThreshold, std::vector<RelationEdge>& vecRelationEdge, Graph& g);

void reduceTransitEdges(const std::vector<Primitive*>& vecPrimitive, std::vector<RelationEdge>& vecRelationEdge, RelationEdge::RelationEdgeType relationType, Graph& g);

#endif // RelationGraph_H
