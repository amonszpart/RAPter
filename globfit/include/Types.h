#ifndef Types_H
#define Types_H

#include <osg/Vec3>
#include <CGAL/Cartesian.h>

typedef CGAL::Cartesian<double>     Kernel;
typedef Kernel::Point_3             Point;
typedef Kernel::Vector_3            Vector;
typedef Kernel::Iso_cuboid_3        Iso_cuboid;

struct RichPoint{
  Point   point;
  Vector  normal;
  double  confidence;
};

template <class Vec3>
class Vec3Caster {
public:
  Vec3Caster(double x, double y, double z):_x(x), _y(y), _z(z) {}
  Vec3Caster(const osg::Vec3& vec3) : _x(vec3.x()), _y(vec3.y()), _z(vec3.z()) {}
  Vec3Caster(const Vec3& vec3) : _x(vec3.x()), _y(vec3.y()), _z(vec3.z())  {}

  operator osg::Vec3() const {return osg::Vec3(_x, _y, _z);}
  operator Vec3() const {return Vec3(_x, _y, _z);}

private:
  double _x, _y, _z;
};

#endif // Types_H


