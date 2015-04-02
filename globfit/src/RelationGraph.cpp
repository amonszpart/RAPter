#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/connected_components.hpp>

#include "Primitive.h"
#include "GlobFit.h"
#include "RelationGraph.h"

class PathDetectedException
{
public:
protected:
private:
};

class path_recorder:public default_bfs_visitor {
public:
    path_recorder(std::vector<GraphVertex>& vecPredecessor, GraphVertex startVertex, GraphVertex endVertex)
        :_vecPredecessor(vecPredecessor), _startVertex(startVertex), _endVertex(endVertex)
    {}

    template < typename Vertex, typename Graph >
    void initialize_vertex(Vertex u, const Graph & g) const
    {
        _vecPredecessor[u] = u;
    }

    template < typename Edge, typename Graph >
    void tree_edge(Edge e, const Graph & g) const
    {
        GraphVertex v = target(e, g);
        if (v == _startVertex) {
            return;
        }
        _vecPredecessor[v] = source(e, g);
        if (v == _endVertex) {
            throw PathDetectedException();
        }
    }
protected:
private:
    std::vector<GraphVertex>& _vecPredecessor;
    GraphVertex _startVertex;
    GraphVertex _endVertex;
};

static void removeIsolatedEdges(Graph& g)
{
    bool bEdgeRemoved;
    do {
        bEdgeRemoved = false;
        for (GraphVertex i = 0, iEnd = num_vertices(g); i < iEnd; ++ i) {
            if (out_degree(i, g) == 1) {
                std::pair<GraphOutEdgeItr, GraphOutEdgeItr> edgeItr = out_edges(i, g);
                GraphVertex v = target(*(edgeItr.first), g);
                if (out_degree(v, g) != 1) {
                    clear_vertex(i, g);
                    bEdgeRemoved = true;
                }
            }
        }
    } while (bEdgeRemoved);
}

static void detectArticulationPoints(std::vector<RelationEdge>& vecRelationEdge, Graph& g)
{
    for (size_t i = 0, iEnd = vecRelationEdge.size(); i < iEnd; ++ i) {
        RelationEdge& relationEdge = vecRelationEdge[i];

        GraphVertex u = relationEdge.getSource().getIdx();
        GraphVertex v = relationEdge.getTarget().getIdx();

        add_edge(u, v, EdgeProp(relationEdge.getType(), relationEdge.getScore()), g);
    }

    removeIsolatedEdges(g);

    std::vector<GraphVertex> vecArticulationPoint;
    articulation_points(g, std::back_inserter(vecArticulationPoint));
    while (!vecArticulationPoint.empty()) {
        GraphVertex nWeakestPoint = 0;
        double nMinWeight = std::numeric_limits<double>::max();
        for (GraphVertex i = 0; i < vecArticulationPoint.size(); ++ i) {
            double nWeight = 0;
            std::pair<GraphOutEdgeItr, GraphOutEdgeItr> edgeItr = out_edges(vecArticulationPoint[i], g);
            for (GraphOutEdgeItr it = edgeItr.first; it != edgeItr.second; it ++) {
                nWeight += g[*it].score;
            }
            if (nWeight < nMinWeight) {
                nMinWeight = nWeight;
                nWeakestPoint = vecArticulationPoint[i];
            }
        }
        std::map<GraphVertex, EdgeProp> mapVertex2Edge;
        std::pair<GraphOutEdgeItr, GraphOutEdgeItr> edgeItr = out_edges(nWeakestPoint, g);
        for (GraphOutEdgeItr it = edgeItr.first; it != edgeItr.second; it ++) {
            mapVertex2Edge[target(*it, g)] = g[*it];
        }
        clear_vertex(nWeakestPoint, g);
        std::vector<GraphVertex> component(num_vertices(g));
        size_t nComponentNum = connected_components(g, &component[0]);
        std::vector<double> vecWeight(nComponentNum, 0.0);
        std::vector<int> vecCount(nComponentNum, 0);
        for (std::map<GraphVertex, EdgeProp>::iterator it = mapVertex2Edge.begin(); it != mapVertex2Edge.end(); it ++) {
            vecWeight[component[it->first]] += it->second.score;
            vecCount[component[it->first]] ++;
        }
        for (size_t i = 0; i < nComponentNum; ++ i) {
            if (vecCount[i] != 0) {
                vecWeight[i] /= vecCount[i];
            }
        }
        size_t nStrongestComponent = std::distance(vecWeight.begin(), std::max_element(vecWeight.begin(), vecWeight.end()));
        for (std::map<GraphVertex, EdgeProp>::iterator it = mapVertex2Edge.begin(); it != mapVertex2Edge.end(); it ++) {
            GraphVertex v = it->first;
            if (component[v] == nStrongestComponent) {
                add_edge(nWeakestPoint, v, mapVertex2Edge[v], g);
            }
        }

        removeIsolatedEdges(g);

        vecArticulationPoint.clear();
        articulation_points(g, std::back_inserter(vecArticulationPoint));
    }

    return;
}

static void filterEdgesToArticulationPoints(std::vector<RelationEdge>& vecRelationEdge, Graph& g)
{
    std::vector<RelationEdge> vecRemainingEdge;
    for (size_t i = 0, iEnd = vecRelationEdge.size(); i < iEnd; ++ i) {
        RelationEdge& relationEdge = vecRelationEdge[i];

        GraphVertex u = relationEdge.getSource().getIdx();
        GraphVertex v = relationEdge.getTarget().getIdx();

        if (edge(u, v, g).second) {
            vecRemainingEdge.push_back(relationEdge);
        }
    }
    vecRelationEdge = vecRemainingEdge;
    return;
}

void biconnectDecompose(std::vector<RelationEdge>& vecRelationEdge, Graph& g)
{
    detectArticulationPoints(vecRelationEdge, g);
    filterEdgesToArticulationPoints(vecRelationEdge, g);

    return;
}

static void filterTransitConstraint(std::vector<RelationEdge>& vecRelationEdge, RelationEdge::RelationEdgeType rt, Graph& g)
{
    std::vector<RelationEdge> vecRemainingEdge;
    for (size_t i = 0, iEnd = vecRelationEdge.size(); i < iEnd; ++ i) {
        RelationEdge& relationEdge = vecRelationEdge[i];

        GraphVertex u = relationEdge.getSource().getIdx();
        GraphVertex v = relationEdge.getTarget().getIdx();

        if (relationEdge.getType() == rt) {
            if (g[u].getParent() == g[v].getParent()) {
                continue;
            }
        }

        vecRemainingEdge.push_back(relationEdge);
    }
    vecRelationEdge = vecRemainingEdge;
    std::sort(vecRelationEdge.begin(), vecRelationEdge.end());

    return;
}

static void collectTransitConstraint(std::vector<RelationEdge>& vecRelationEdge, RelationEdge::RelationEdgeType rt, Graph& g)
{
    for (GraphVertex i = 0, iEnd = num_vertices(g); i < iEnd; ++ i) {
        if (g[i].getParent() != i) {
            const RelationVertex& targetVertex = g[g[i].getParent()];
            const RelationVertex& sourceVertex = g[i];
            vecRelationEdge.push_back(RelationEdge(rt, sourceVertex, targetVertex));
        }
    }

    return;
}

static void filterOrientationConstraint(std::vector<RelationEdge>& vecRelationEdge, Graph& g)
{
    filterTransitConstraint(vecRelationEdge, RelationEdge::RET_PARALLEL, g);

    std::vector<RelationEdge> vecRemainingEdge;
    for (size_t i = 0, iEnd = vecRelationEdge.size(); i < iEnd; ++ i) {
        RelationEdge& relationEdge = vecRelationEdge[i];

        GraphVertex u = relationEdge.getSource().getIdx();
        GraphVertex v = relationEdge.getTarget().getIdx();

        if (relationEdge.getType() == RelationEdge::RET_ORTHOGONAL) {
            if (edge(g[u].getParent(), g[v].getParent(), g).second) {
                continue;
            }
        }
        vecRemainingEdge.push_back(relationEdge);
    }
    vecRelationEdge = vecRemainingEdge;
    std::sort(vecRelationEdge.begin(), vecRelationEdge.end());

    return;
}

static void collectOrientationConstraint(std::vector<RelationEdge>& vecRelationEdge, Graph& g)
{
    collectTransitConstraint(vecRelationEdge, RelationEdge::RET_PARALLEL, g);

    std::pair<GraphEdgeItr, GraphEdgeItr> edgeItr = edges(g);

    for (GraphEdgeItr it = edgeItr.first; it != edgeItr.second; it ++) {
        GraphVertex u = source(*it, g);
        GraphVertex v = target(*it, g);

        const RelationVertex& sourceVertex = g[u];
        const RelationVertex& targetVertex = g[v];

        vecRelationEdge.push_back(RelationEdge(RelationEdge::RET_ORTHOGONAL, sourceVertex, targetVertex));
    }

    return;
}

static bool findBridge(GraphVertex u, GraphVertex v, const Graph& g, GraphVertex &bridge)
{
	std::set<GraphVertex> setBridge;
	std::pair<GraphOutEdgeItr, GraphOutEdgeItr> edgeItr = out_edges(u, g);
	for (GraphOutEdgeItr it = edgeItr.first; it != edgeItr.second; it ++) {
		setBridge.insert(target(*it, g));
	}

	edgeItr = out_edges(v, g);
	for (GraphOutEdgeItr it = edgeItr.first; it != edgeItr.second; it ++) {
		bridge = target(*it, g);
		if (setBridge.find(bridge) != setBridge.end()) {
			return true;
		}
	}

	return false;
}

static std::pair<GraphVertex, GraphVertex> findParalell(const std::vector<GraphVertex>& vecLoop, const std::vector<Primitive*>& vecPrimitive, bool bForce, double angleThreshold)
{
	std::pair<GraphVertex, GraphVertex> para;
	double nMinAngle = std::numeric_limits<double>::max();

	int nLoopSize = vecLoop.size();
	for (int i = 0; i < nLoopSize; ++ i) {
		Vector sNormal;
		if (!vecPrimitive[vecLoop[i]]->getNormal(sNormal)) {
			continue;
		}
		for (int j = i + 2; j < nLoopSize; ++ j) {
			if (i == 0 && j == nLoopSize-1) {
				continue;
			}

			Vector tNormal;
			if (!vecPrimitive[vecLoop[j]]->getNormal(tNormal)) {
				continue;
			}

			double nAngle = std::acos(std::abs(sNormal*tNormal));
			if (nAngle < nMinAngle) {
				nMinAngle = nAngle;
				para.first = vecLoop[i];
				para.second = vecLoop[j];
			}
		}
	}

	if (!bForce && nMinAngle > angleThreshold) {
		return std::pair<GraphVertex, GraphVertex>(0, 0);
	}

	return para;
}

static bool willProduceLoop(std::vector<GraphVertex>& vecPredecessor, GraphVertex v, GraphVertex u, Graph& g)
{
    path_recorder pathRecorder(vecPredecessor, v, u);
    try {
        breadth_first_search(g, v, visitor(pathRecorder));
        // no path, so after we add edge(u, v), there will be no loop
        return false;
    } catch (PathDetectedException& e) {
        // there is a path between u and v, and we added edge(u, v), so there will be loop
        // and edge(u,v) is redundant
        return true;
    }
}

static void addParaEdge(GraphVertex u, GraphVertex v,
    const std::vector<Primitive*>& vecPrimitive,
    std::vector<GraphVertex>& vecPredecessor,
    double angleThreshold,
    Graph& g);

static void addOrthEdge(GraphVertex u, GraphVertex v,
    const std::vector<Primitive*>& vecPrimitive,
    std::vector<GraphVertex>& vecPredecessor,
    double angleThreshold,
    Graph& g)
{
	u = g[u].getParent();
	v = g[v].getParent();
	// check if edge exists
	if (edge(u, v, g).second) {
		return;
	}

    if (!willProduceLoop(vecPredecessor, v, u, g)) {
        add_edge(u, v, EdgeProp(RelationEdge::RET_ORTHOGONAL), g);
    }
    else {
        add_edge(u, v, EdgeProp(RelationEdge::RET_ORTHOGONAL), g);

        // there is a path between u and v, and we added edge(u, v), so there will be loop
        std::vector<GraphVertex> vecLoop;
        GraphVertex current = u;
        vecLoop.push_back(current);
        while (vecPredecessor[current] != current) {
            current = vecPredecessor[current];
            vecLoop.push_back(current);
        }

        int nLoopSize = vecLoop.size();
        if (nLoopSize == 3) {
            GraphVertex opposite2 = 0;
            GraphVertex opposite0 = 0;

            bool bShare01 = findBridge(vecLoop[0], vecLoop[1], g, opposite2);
            bool bShare12 = findBridge(vecLoop[1], vecLoop[2], g, opposite0);

            if (bShare01) {
                addParaEdge(opposite2, vecLoop[2], vecPrimitive, vecPredecessor, angleThreshold, g);
            }
            if (bShare12) {
                addParaEdge(opposite0, vecLoop[0], vecPrimitive, vecPredecessor, angleThreshold, g);
            }
        } else {
            std::pair<GraphVertex, GraphVertex> para = findParalell(vecLoop, vecPrimitive, (nLoopSize == 4), angleThreshold);
            if (para.first != para.second) {
                addParaEdge(para.first, para.second, vecPrimitive, vecPredecessor, angleThreshold, g);
            }
        }
    }

	return;
}

static void addParaEdge(GraphVertex u, GraphVertex v,
    const std::vector<Primitive*>& vecPrimitive,
    std::vector<GraphVertex>& vecPredecessor,
    double angleThreshold,
    Graph& g)
{
	u = g[u].getParent();
	v = g[v].getParent();

	if (u == v) {
		return;
	}

	if (degree(u, g) != 0) {
		std::swap(u, v);
	}

	// path compression
	GraphVertex nVertexNum = num_vertices(g);
	std::vector<GraphVertex> vecStar;
	for (GraphVertex i = 0; i < nVertexNum; ++ i) {
		if (g[i].getParent() == u) {
			g[i].setParent(v);
			vecStar.push_back(i);
		}
	}

	for (std::vector<GraphVertex>::iterator it = vecStar.begin(); it != vecStar.end(); it ++) {
		GraphVertex i = *it;
		std::vector<GraphVertex> vecBridge;
		std::pair<GraphOutEdgeItr, GraphOutEdgeItr> edgeItr = out_edges(i, g);
		for (GraphOutEdgeItr it = edgeItr.first; it != edgeItr.second; it ++) {
			vecBridge.push_back(target(*it, g));
		}

		int nBridgeNum = vecBridge.size();
		for (int j = 0; j < nBridgeNum; ++ j) {
			remove_edge(i, vecBridge[j], g);
		}
		for (int j = 0; j < nBridgeNum; ++ j) {
			addOrthEdge(vecBridge[j], v, vecPrimitive, vecPredecessor, angleThreshold, g);
		}
	}

	return;
}

void reduceParaOrthEdges(const std::vector<Primitive*>& vecPrimitive, double angleThreshold, std::vector<RelationEdge>& vecRelationEdge, Graph& g)
{
    filterOrientationConstraint(vecRelationEdge, g);

    if (vecRelationEdge.size() == 0) {
        collectOrientationConstraint(vecRelationEdge, g);
        return;
    }

    while (!vecRelationEdge.empty()) {
        RelationEdge& relationEdge = vecRelationEdge.back();
        vecRelationEdge.pop_back();
        GraphVertex u = relationEdge.getSource().getIdx();
        GraphVertex v = relationEdge.getTarget().getIdx();
        std::vector<GraphVertex> vecPredecessor(vecPrimitive.size());
        if (relationEdge.getType() == RelationEdge::RET_PARALLEL) {
            addParaEdge(u, v, vecPrimitive, vecPredecessor, angleThreshold, g);
        } else if (relationEdge.getType() == RelationEdge::RET_ORTHOGONAL) {
            addOrthEdge(u, v, vecPrimitive, vecPredecessor, angleThreshold, g);
        }

        filterOrientationConstraint(vecRelationEdge, g);
    }

    collectOrientationConstraint(vecRelationEdge, g);

    return;
}

void reduceTransitEdges(const std::vector<Primitive*>& vecPrimitive, std::vector<RelationEdge>& vecRelationEdge, RelationEdge::RelationEdgeType relationType, Graph& g)
{
    std::sort(vecRelationEdge.begin(), vecRelationEdge.end());

    std::map<GraphVertex, RelationVertex> mapVertex;
    for (size_t i = 0, iEnd = vecRelationEdge.size(); i < iEnd; ++ i) {
        GraphVertex u = vecRelationEdge[i].getTarget().getIdx();
        GraphVertex v = vecRelationEdge[i].getSource().getIdx();

        mapVertex[u] = vecRelationEdge[i].getTarget();
        mapVertex[v] = vecRelationEdge[i].getSource();

        add_edge(u, v, EdgeProp(relationType), g);
    }

    std::vector<GraphVertex> component(num_vertices(g));
    size_t numComponent = connected_components(g, &component[0]);

    std::map<size_t, std::vector<size_t> > mapComponent;
    for (size_t i = 0, iEnd = component.size(); i < iEnd; ++ i) {
        mapComponent[component[i]].push_back(i);
    }

    vecRelationEdge.clear();
    for (std::map<size_t, std::vector<size_t> >::const_iterator it = mapComponent.begin();
        it != mapComponent.end();
        ++ it) {
        const std::vector<size_t>& vecComponent = it->second;
        if (vecComponent.size() < 2) {
            continue;
        }
        GraphVertex u = vecComponent[0];
        for (size_t i = 1, iEnd = vecComponent.size(); i < iEnd; ++ i) {
            GraphVertex v = vecComponent[i];
            RelationEdge relationEdge(relationType, mapVertex[v], mapVertex[u], 1.0);
            GlobFit::computeEdgeScore(relationEdge, vecPrimitive);
            vecRelationEdge.push_back(relationEdge);
        }
    }
    std::sort(vecRelationEdge.begin(), vecRelationEdge.end());

    if (vecRelationEdge.size() == 0) {
        return;
    }

    if (vecRelationEdge[0].getType() != RelationEdge::RET_EQUAL_ANGLE && vecRelationEdge[0].getType() != RelationEdge::RET_EQUAL_LENGTH) {
        return;
    }

    g.clear();
    for (size_t i = 0, iEnd = vecPrimitive.size(); i < iEnd; ++ i) {
        RelationVertex relationVertex(i, vecPrimitive[i]->getIdx());
        relationVertex.setParent(i);
        add_vertex(relationVertex, g);
    }

    std::vector<RelationEdge> vecRemainingEdge;
    std::vector<GraphVertex> vecPredecessor(vecPrimitive.size());
    for (size_t i = 0, iEnd = vecRelationEdge.size(); i < iEnd; ++ i) {
        const RelationVertex& source = vecRelationEdge[i].getSource();
        const RelationVertex& target = vecRelationEdge[i].getTarget();

        bool sourceProduceLoop = willProduceLoop(vecPredecessor, source.getPrimitiveIdx1(), source.getPrimitiveIdx2(), g);
        bool targetProduceLoop = willProduceLoop(vecPredecessor, target.getPrimitiveIdx1(), target.getPrimitiveIdx2(), g);

        // this is too strong for avoiding conflicts, some compatible cases may be removed
        // TODO: deduce the right rules for edges with 4 primitives involved
        if (sourceProduceLoop || targetProduceLoop) {
            continue;
        }

        add_edge(source.getPrimitiveIdx1(), source.getPrimitiveIdx2(), g);
        add_edge(target.getPrimitiveIdx1(), target.getPrimitiveIdx2(), g);

        vecRemainingEdge.push_back(vecRelationEdge[i]);
    }

    vecRelationEdge = vecRemainingEdge;

    return;
}