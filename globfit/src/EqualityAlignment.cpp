#include "Types.h"
#include "Primitive.h"
#include "RelationGraph.h"

#include "GlobFit.h"

bool GlobFit::equalLengthAlignment(double equalLengthThreshold)
{
  _vecCoplanarCollapse.resize(_vecPrimitive.size());
  for (size_t i = 0, iEnd = _vecCoplanarCollapse.size(); i < iEnd; ++ i) {
    _vecCoplanarCollapse[i] = i;
  }

  for (size_t i = 0, iEnd = _vecDistanceEdge.size(); i < iEnd; ++ i) {
    const RelationEdge& relationEdge = _vecDistanceEdge[i];
    if (relationEdge.getType() == RelationEdge::RET_COPLANAR) {
      _vecCoplanarCollapse[relationEdge.getSource().getPrimitiveIdx1()] = relationEdge.getTarget().getPrimitiveIdx1();
    }
  }

  std::vector<std::set<size_t> > vecParallelPlaneGroup(_vecPrimitive.size());

  for (size_t i = 0, iEnd = _vecParallelCollapse.size(); i < iEnd; ++ i) {
    if (_vecPrimitive[i]->getType() == Primitive::PT_PLANE) {
      vecParallelPlaneGroup[_vecParallelCollapse[i]].insert(_vecCoplanarCollapse[i]);
    }
  }

  Graph lengthGraph;
  std::vector<RelationVertex> vecLengthVertex;
  std::vector<double> vecLength;
  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    std::vector<size_t> vecParallelPlane(vecParallelPlaneGroup[i].begin(), vecParallelPlaneGroup[i].end());
    for (size_t j = 0, jEnd = vecParallelPlane.size(); j < jEnd; ++ j) {
      for (size_t k = j+1, kEnd = vecParallelPlane.size(); k < kEnd; ++ k) {
        size_t idx1 = vecParallelPlane[j];
        size_t idx2 = vecParallelPlane[k];

        RelationVertex relationVertex(vecLengthVertex.size(), idx1, idx2);
        relationVertex.setParent(vecLengthVertex.size());

        vecLengthVertex.push_back(relationVertex);
        add_vertex(relationVertex, lengthGraph);

        double d1, d2;
        _vecPrimitive[idx1]->getDistance(d1);
        _vecPrimitive[idx2]->getDistance(d2);

        Vector normal1, normal2;
        _vecPrimitive[idx1]->getNormal(normal1);
        _vecPrimitive[idx2]->getNormal(normal2);

        vecLength.push_back(computeLength(normal1, normal2, d1, d2));
      }
    }
  }

  std::vector<RelationEdge> vecLengthEdge;
  for (size_t i = 0, iEnd = vecLengthVertex.size(); i < iEnd; ++ i) {
    for (size_t j = i+1, jEnd = vecLengthVertex.size(); j < jEnd; ++ j) {
      double lengthScore = computeEqualLengthScore(vecLength[i], vecLength[j]);
      if(lengthScore > equalLengthThreshold) {
        vecLengthEdge.push_back(RelationEdge(RelationEdge::RET_EQUAL_LENGTH, vecLengthVertex[i], vecLengthVertex[j], lengthScore));
      }
    }
  }

  reduceTransitEdges(_vecPrimitive, vecLengthEdge, RelationEdge::RET_EQUAL_LENGTH, lengthGraph);

  if (vecLengthEdge.size() == 0) {
    return true;
  }

  _vecDistanceEdge.insert(_vecDistanceEdge.end(), vecLengthEdge.begin(), vecLengthEdge.end());
  while(!solve(_vecDistanceEdge, RelationEdge::RET_EQUAL_LENGTH, "EqualLength")) {
    _vecDistanceEdge.pop_back();
  }

  return true;
}

bool GlobFit::equalRadiusAlignment(double equalRadiusThreshold)
{
  Graph equalRadiusGraph;
  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    RelationVertex relationVertex(i, _vecPrimitive[i]->getIdx());
    relationVertex.setParent(i);
    add_vertex(relationVertex, equalRadiusGraph);
  }

  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    Primitive* p1 = _vecPrimitive[i];
    double radius1;
    if (!p1->getRadius(radius1)) {
      continue;
    }

    for (size_t j = i + 1, jEnd = iEnd; j < jEnd; ++ j) {
      Primitive* p2 = _vecPrimitive[j];
      double radius2;
      if (!p2->getRadius(radius2)) {
        continue;
      }

      double equalRadiusScore = computeEqualRadiusScore(radius1, radius2);
      if (equalRadiusScore > equalRadiusThreshold) {
        _vecRadiusEdge.push_back(RelationEdge(RelationEdge::RET_EQUAL_RADIUS, equalRadiusGraph[i], equalRadiusGraph[j], equalRadiusScore));
      }
    }
  }

  reduceTransitEdges(_vecPrimitive, _vecRadiusEdge, RelationEdge::RET_EQUAL_RADIUS, equalRadiusGraph);

  return solve(_vecRadiusEdge, RelationEdge::RET_EQUAL_RADIUS, "EqualRadius");
}

bool GlobFit::equalityAlignment(double equalLengthThreshld, double equalRadiusThreshold)
{
  equalLengthThreshld = equalLengthThreshld > 0 ? -equalLengthThreshld:equalLengthThreshld;
  equalRadiusThreshold = equalRadiusThreshold > 0 ? -equalRadiusThreshold:equalRadiusThreshold;

  if(!equalLengthAlignment(equalLengthThreshld)) {
    return false;
  }

  _vecRadiusEdge.clear();
  if(!equalRadiusAlignment(equalRadiusThreshold)) {
    return false;
  }

  return true;
}