#ifndef GlobFit_H
#define GlobFit_H

#include <string>
#include <vector>

#include "RelationEdge.h"
#include "Types.h"
#include "CoreExports.h"

struct  RichPoint;
class   Primitive;

namespace osg{
  class Node;
}

class CORE_EXPORTS GlobFit
{
public:
  GlobFit(void);
  ~GlobFit(void);

  bool load(const std::string& filename);
  bool save(const std::string& filename) const;

  bool createMatlabArraies();
  void destoryMatlabArraies();

  bool orientationAlignment(double paraOrthThreshold, double equalAngleThreshold);
  bool placementAlignment(double coaxialThreshold, double coplanarThreshold);
  bool equalityAlignment(double equalLengthThreshold, double equalRadiusThreshold);

  const std::vector<RichPoint*>&  getPointSet() const {return _vecPointSet;}
  const std::vector<Primitive*>&  getPrimitive() const {return _vecPrimitive;}

  std::pair<osg::Node*, osg::Node*> convertPointsToGeometry() const;
  osg::Node* convertPrimitivesToGeometry(const std::string& title) const;

  static double computeAngle(const Vector& normal1, const Vector& normal2);
  static double computeLength(const Vector& normal1, const Vector& normal2, double d1, double d2);
  static double computeParallelScore(const Vector& normal1, const Vector& normal2);
  static double computeOrthogonalScore(const Vector& normal1, const Vector& normal2);
  static double computeEqualAngleScore(double angle1, double angle2);
  static double computeCoaxialScore(Vector normal1, Point point1, Point point2);
  static double computeCoplanarScore(double d1, double d2);
  static double computeEqualRadiusScore(double r1, double r2);
  static double computeEqualLengthScore(double l1, double l2);

  static void computeEdgeScore(RelationEdge& relationEdge, const std::vector<Primitive*>& vecPrimitive);

protected:
  bool solve(std::vector<RelationEdge>& vecEdge, RelationEdge::RelationEdgeType currentStage, const std::string& stageName);
  void dumpData(const std::vector<RelationEdge>& vecEdge, const std::string& stageName);

  bool paraOrthAlignment(double orientationThreshold);
  bool equalAngleAlignment(double orientationThreshold);

  bool coaxialAlignment(double coaxialThreshold);
  bool coplanarAlignment(double coplanarThreshold);

  bool equalLengthAlignment(double equalLengthThreshold);
  bool equalRadiusAlignment(double equalRadiusThreshold);

private:
  std::vector<RichPoint*>     _vecPointSet;
  std::vector<Primitive*>     _vecPrimitive;

  std::vector<RelationEdge>   _vecNormalEdge;
  std::vector<RelationEdge>   _vecPointEdge;
  std::vector<RelationEdge>   _vecDistanceEdge;
  std::vector<RelationEdge>   _vecRadiusEdge;
  std::vector<size_t>         _vecParallelCollapse;
  std::vector<size_t>         _vecCoplanarCollapse;
};

#endif // GlobFit_H