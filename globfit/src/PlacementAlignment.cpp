#include "Types.h"
#include "Primitive.h"
#include "RelationGraph.h"

#include "GlobFit.h"

bool GlobFit::coaxialAlignment(double coaxialThreshold)
{
  Graph coaxialGraph;
  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    RelationVertex relationVertex(i, _vecPrimitive[i]->getIdx());
    relationVertex.setParent(i);
    add_vertex(relationVertex, coaxialGraph);
  }

  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    Primitive* p1 = _vecPrimitive[i];
    Vector normal1;
    if (!p1->getNormal(normal1)) {
      continue;
    }
    Point point1;
    if (!p1->getCenter(point1)) {
      continue;
    }

    for (size_t j = i + 1, jEnd = iEnd; j < jEnd; ++ j) {
      Primitive* p2 = _vecPrimitive[j];
      Vector normal2;
      if (!p2->getNormal(normal2)) {
        continue;
      }
      Point point2;
      if (!p2->getCenter(point2)) {
        continue;
      }
      if(_vecParallelCollapse[i] != _vecParallelCollapse[j]) {
        continue;
      }

      double coaxialScore = computeCoaxialScore(normal1, point1, point2);
      if (coaxialScore > coaxialThreshold) {
        _vecPointEdge.push_back(RelationEdge(RelationEdge::RET_COAXIAL, coaxialGraph[i], coaxialGraph[j], coaxialScore));
      }
    }
  }

  reduceTransitEdges(_vecPrimitive, _vecPointEdge, RelationEdge::RET_COAXIAL, coaxialGraph);

  std::vector<size_t> vecCoaxialCollapse(_vecPrimitive.size());
  for (size_t i = 0, iEnd = vecCoaxialCollapse.size(); i < iEnd; ++ i) {
    vecCoaxialCollapse[i] = i;
  }
  for (size_t i = 0, iEnd = _vecPointEdge.size(); i < iEnd; ++ i) {
    const RelationEdge& relationEdge = _vecPointEdge[i];
    vecCoaxialCollapse[relationEdge.getSource().getPrimitiveIdx1()] = relationEdge.getTarget().getPrimitiveIdx1();
  }

  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    Primitive* p1 = _vecPrimitive[i];
    Vector normal1;
    if (p1->getNormal(normal1)) {
      continue;
    }
    Point point1;
    if (!p1->getCenter(point1)) {
      continue;
    }

    std::map<size_t, double> mapTarget;
    for (size_t j = i + 1, jEnd = iEnd; j < jEnd; ++ j) {
      Primitive* p2 = _vecPrimitive[j];
      Vector normal2;
      if (!p2->getNormal(normal2)) {
        continue;
      }
      Point point2;
      if (!p2->getCenter(point2)) {
        continue;
      }

      double coaxialScore = computeCoaxialScore(normal2, point1, point2);
      if (coaxialScore > coaxialThreshold) {
        std::map<size_t, double>::iterator it = mapTarget.find(vecCoaxialCollapse[j]);
        if (it == mapTarget.end()) {
          mapTarget[vecCoaxialCollapse[j]] = coaxialScore;
        } else if(it->second < coaxialScore) {
          it->second = coaxialScore;
        }
      }
    }
    for (std::map<size_t, double>::iterator it = mapTarget.begin(); it != mapTarget.end(); ++ it) {
      _vecPointEdge.push_back(RelationEdge(RelationEdge::RET_COAXIAL, coaxialGraph[i], coaxialGraph[it->first], it->second));
    }
  }

  return solve(_vecPointEdge, RelationEdge::RET_COAXIAL, "CoAxial");
}

bool GlobFit::coplanarAlignment(double coplanarThreshold)
{
  Graph coplanarGraph;
  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    RelationVertex relationVertex(i, _vecPrimitive[i]->getIdx());
    relationVertex.setParent(i);
    add_vertex(relationVertex, coplanarGraph);
  }

  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    Primitive* p1 = _vecPrimitive[i];
    if (p1->getType() != Primitive::PT_PLANE) {
      continue;
    }

    double d1;
    p1->getDistance(d1);

    for (size_t j = i + 1, jEnd = iEnd; j < jEnd; ++ j) {
      Primitive* p2 = _vecPrimitive[j];
      Primitive* p1 = _vecPrimitive[i];
      if (p2->getType() != Primitive::PT_PLANE) {
        continue;
      }

      if (_vecParallelCollapse[i] != _vecParallelCollapse[j]) {
        continue;
      }

      double d2;
      p2->getDistance(d2);
      double coplanarScore = computeCoplanarScore(d1, d2);
      if (coplanarScore > coplanarThreshold) {
        _vecDistanceEdge.push_back(RelationEdge(RelationEdge::RET_COPLANAR, coplanarGraph[i], coplanarGraph[j], coplanarScore));
      }
    }
  }

  reduceTransitEdges(_vecPrimitive, _vecDistanceEdge, RelationEdge::RET_COPLANAR, coplanarGraph);

  return solve(_vecDistanceEdge, RelationEdge::RET_COPLANAR, "CoPlanar");
}

bool GlobFit::placementAlignment(double coaxialThreshold, double coplanarThreshold)
{
  coaxialThreshold = coaxialThreshold > 0 ? -coaxialThreshold:coaxialThreshold;
  coplanarThreshold = coplanarThreshold > 0 ? -coplanarThreshold:coplanarThreshold;

  _vecPointEdge.clear();
  if(!coaxialAlignment(coaxialThreshold)) {
    return false;
  }

  _vecDistanceEdge.clear();
  if(!coplanarAlignment(coplanarThreshold)) {
    return false;
  }

  return true;
}