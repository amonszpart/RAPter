#define _USE_MATH_DEFINES
#include <cmath>

#include "Types.h"
#include "Primitive.h"
#include "RelationGraph.h"

#include "GlobFit.h"

bool GlobFit::paraOrthAlignment(double paraOrthThreshold)
{
  Graph paraOrthGraph;
  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    RelationVertex relationVertex(i, _vecPrimitive[i]->getIdx());
    relationVertex.setParent(i);
    add_vertex(relationVertex, paraOrthGraph);
  }

  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    Primitive* p1 = _vecPrimitive[i];
    Vector normal1;
    if (!p1->getNormal(normal1)) {
      continue;
    }

    for (size_t j = i + 1, jEnd = iEnd; j < jEnd; ++ j) {
      Primitive* p2 = _vecPrimitive[j];
      Vector normal2;
      if (!p2->getNormal(normal2)) {
        continue;
      }

      double parallelScore = computeParallelScore(normal1, normal2);
      if (parallelScore > paraOrthThreshold) {
        _vecNormalEdge.push_back(RelationEdge(RelationEdge::RET_PARALLEL, paraOrthGraph[i], paraOrthGraph[j], parallelScore));
        continue;
      }

      double orthogonalScore = computeOrthogonalScore(normal1, normal2);
      if (orthogonalScore > paraOrthThreshold) {
        _vecNormalEdge.push_back(RelationEdge(RelationEdge::RET_ORTHOGONAL, paraOrthGraph[i], paraOrthGraph[j], orthogonalScore));
        continue;
      }
    }
  }

  Graph biconnectGraph;
  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    add_vertex(biconnectGraph);
  }
  biconnectDecompose(_vecNormalEdge, biconnectGraph);

  reduceParaOrthEdges(_vecPrimitive, paraOrthThreshold, _vecNormalEdge, paraOrthGraph);

  return solve(_vecNormalEdge, RelationEdge::RET_ORTHOGONAL, "ParaOrth");
}

bool GlobFit::equalAngleAlignment(double equalAngleThreshold)
{
  _vecParallelCollapse.resize(_vecPrimitive.size());
  for (size_t i = 0, iEnd = _vecParallelCollapse.size(); i < iEnd; ++ i) {
    _vecParallelCollapse[i] = i;
  }

  std::set<RelationVertex> setOrthVertex;
  for (size_t i = 0, iEnd = _vecNormalEdge.size(); i < iEnd; ++ i) {
    const RelationEdge& relationEdge = _vecNormalEdge[i];
    if (relationEdge.getType() == RelationEdge::RET_PARALLEL) {
      _vecParallelCollapse[relationEdge.getSource().getPrimitiveIdx1()] = relationEdge.getTarget().getPrimitiveIdx1();
    } else {
      size_t sourceIdx = relationEdge.getSource().getPrimitiveIdx1();
      size_t targetIdx = relationEdge.getTarget().getPrimitiveIdx1();
      if (sourceIdx > targetIdx) {
        std::swap(sourceIdx, targetIdx);
      }
      setOrthVertex.insert(RelationVertex(0, sourceIdx, targetIdx));
    }
  }

  std::vector<size_t> vecRoot;
  for (size_t i = 0, iEnd = _vecParallelCollapse.size(); i < iEnd; ++ i) {
    Vector normal;
    if (_vecParallelCollapse[i] == i && _vecPrimitive[i]->getNormal(normal)) {
      vecRoot.push_back(i);
    }
  }

  Graph angleGraph;
  std::vector<RelationVertex> vecAngleVertex;
  std::vector<double> vecAngle;
  for (size_t i = 0, iEnd = vecRoot.size(); i < iEnd; ++ i) {
    for (size_t j = i+1, jEnd = vecRoot.size(); j < jEnd; ++ j) {
      size_t idx1 = vecRoot[i];
      size_t idx2 = vecRoot[j];
      RelationVertex relationVertex(0, idx1, idx2);
      // if it's already orthogonal, ignore it
      if (setOrthVertex.find(relationVertex) != setOrthVertex.end()) {
        continue;
      }

      relationVertex.setIdx(vecAngleVertex.size());
      relationVertex.setParent(vecAngleVertex.size());

      add_vertex(relationVertex, angleGraph);
      vecAngleVertex.push_back(relationVertex);

      Vector sourceNormal, targetNormal;
      _vecPrimitive[idx1]->getNormal(sourceNormal);
      _vecPrimitive[idx2]->getNormal(targetNormal);
      vecAngle.push_back(computeAngle(sourceNormal, targetNormal));
    }
  }

  std::vector<RelationEdge> vecAngleEdge;
  for (size_t i = 0, iEnd = vecAngleVertex.size(); i < iEnd; ++ i) {
    for (size_t j = i+1, jEnd = vecAngleVertex.size(); j < jEnd; ++ j) {
      double angleScore = computeEqualAngleScore(vecAngle[i], vecAngle[j]);
      if(angleScore > equalAngleThreshold) {
        vecAngleEdge.push_back(RelationEdge(RelationEdge::RET_EQUAL_ANGLE, vecAngleVertex[i], vecAngleVertex[j], angleScore));
      }
    }
  }
  reduceTransitEdges(_vecPrimitive, vecAngleEdge, RelationEdge::RET_EQUAL_ANGLE, angleGraph);

  if (vecAngleEdge.size() == 0) {
    return true;
  }

  _vecNormalEdge.insert(_vecNormalEdge.end(), vecAngleEdge.begin(), vecAngleEdge.end());
  while(!solve(_vecNormalEdge, RelationEdge::RET_EQUAL_ANGLE, "EqualAngle")) {
    _vecNormalEdge.pop_back();
  }

  return true;
}

bool GlobFit::orientationAlignment(double paraOrthThreshold, double equalAngleThreshold)
{
  paraOrthThreshold = -std::abs(paraOrthThreshold*M_PI/180);
  equalAngleThreshold = -std::abs(equalAngleThreshold*M_PI/180);

  _vecNormalEdge.clear();

  if(!paraOrthAlignment(paraOrthThreshold)) {
    return false;
  }

  if(!equalAngleAlignment(equalAngleThreshold)) {
    return false;
  }

  return true;
}