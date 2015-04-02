#define _USE_MATH_DEFINES
#include <cmath>

#include "Primitive.h"
#include "GlobFit.h"

double GlobFit::computeAngle(const Vector& normal1, const Vector& normal2)
{
  double cos = std::abs(normal1*normal2)/std::sqrt(normal1.squared_length()*normal2.squared_length());
  cos = std::min(cos, 1.0);
  cos = std::max(cos, 0.0);
  return std::acos(cos);
}

double GlobFit::computeLength(const Vector& normal1, const Vector& normal2, double d1, double d2)
{
  return normal1*normal2 < 0? abs(d1+d2):abs(d1-d2);
}

double GlobFit::computeParallelScore(const Vector& normal1, const Vector& normal2)
{
  return -computeAngle(normal1, normal2);
}

double GlobFit::computeOrthogonalScore(const Vector& normal1, const Vector& normal2)
{
  return computeAngle(normal1, normal2)-M_PI_2;
}

double GlobFit::computeEqualAngleScore(double angle1, double angle2)
{
  return -std::abs(angle1-angle2);
}

double GlobFit::computeCoaxialScore(Vector normal1, Point point1, Point point2)
{
  Kernel::Line_3 line(point1, normal1);

  return -sqrt(CGAL::squared_distance(line, point2));
}

double GlobFit::computeCoplanarScore(double d1, double d2)
{
  return -abs(d1-d2);
}

double GlobFit::computeEqualRadiusScore(double r1, double r2)
{
  return -abs(r1-r2);
}

double GlobFit::computeEqualLengthScore(double l1, double l2)
{
  return -abs(l1-l2);
}

void GlobFit::computeEdgeScore(RelationEdge& relationEdge, const std::vector<Primitive*>& vecPrimitive)
{
  if (relationEdge.getType() == RelationEdge::RET_PARALLEL) {
    Vector normal1, normal2;
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getNormal(normal1);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getNormal(normal2);
    relationEdge.setScore(computeParallelScore(normal1, normal2));
  } else if (relationEdge.getType() == RelationEdge::RET_ORTHOGONAL) {
    Vector normal1, normal2;
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getNormal(normal1);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getNormal(normal2);
    relationEdge.setScore(computeOrthogonalScore(normal1, normal2));
  } else if (relationEdge.getType() == RelationEdge::RET_EQUAL_ANGLE) {
    Vector normal1, normal2, normal3, normal4;
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getNormal(normal1);
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx2()]->getNormal(normal2);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getNormal(normal3);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx2()]->getNormal(normal4);
    relationEdge.setScore(
      computeEqualAngleScore(
      computeAngle(normal1, normal2), computeAngle(normal3, normal4)));
  } else if (relationEdge.getType() == RelationEdge::RET_COAXIAL) {
    Vector normal;
    Point point1, point2;
    if(!vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getNormal(normal)) {
      vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getNormal(normal);
    }
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getCenter(point1);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getCenter(point2);
    relationEdge.setScore(computeCoaxialScore(normal, point1, point2));
  } else if (relationEdge.getType() == RelationEdge::RET_COPLANAR) {
    double d1, d2;
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getDistance(d1);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getDistance(d2);
    relationEdge.setScore(computeCoplanarScore(d1, d2));
  } else if (relationEdge.getType() == RelationEdge::RET_EQUAL_LENGTH) {
    Vector normal1, normal2, normal3, normal4;
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getNormal(normal1);
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx2()]->getNormal(normal2);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getNormal(normal3);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx2()]->getNormal(normal4);

    double d1, d2, d3, d4;
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getDistance(d1);
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx2()]->getDistance(d2);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getDistance(d3);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx2()]->getDistance(d4);
    relationEdge.setScore(
      computeEqualLengthScore(
      computeLength(normal1, normal2, d1, d2), computeLength(normal3, normal4, d3, d4)));
  } else if (relationEdge.getType() == RelationEdge::RET_EQUAL_RADIUS) {
    double r1, r2;
    vecPrimitive[relationEdge.getSource().getPrimitiveIdx1()]->getRadius(r1);
    vecPrimitive[relationEdge.getTarget().getPrimitiveIdx1()]->getRadius(r2);
    relationEdge.setScore(computeEqualRadiusScore(r1, r2));
  }

  return;
}