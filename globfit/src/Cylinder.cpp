#include <sstream>

#include <osg/Geode>
#include <osg/Shape>
#include <osg/ShapeDrawable>

#include "Cylinder.h"


Cylinder::Cylinder(const std::vector<RichPoint*>& vecPointSet) : Primitive(vecPointSet, PT_CYLINDER)
{
}


Cylinder::~Cylinder(void)
{
}

osg::Group* Cylinder::toGeometry(const osg::Vec4& color)
{
    osg::Group* group = new osg::Group;
    if (_vecPointIdx.size() < 2) {
        return group;
    }

    Kernel::Line_3 axis(_point, _normal);

    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    for (size_t i = 0, iEnd = _vecPointIdx.size(); i < iEnd; ++ i) {
        Point projection = axis.projection(_vecPointSet[_vecPointIdx[i]]->point);
        Vector vector(_point, projection);
        double d = vector*_normal;
        max = max<d ? d:max;
        min = min>d ? d:min;
    }

    Point top = _point + max*_normal;
    Point bottom = _point + min*_normal;
    Point center = CGAL::midpoint(top, bottom);
    osg::Vec3 offset = Vec3Caster<Vector>(top-bottom);
    osg::Vec3 zAxis(0.0, 0.0, 1.0);
    osg::Vec3 rotation = zAxis^offset;
    float angle = acos((zAxis*offset)/offset.length());
    osg::Cylinder* cylinder = new osg::Cylinder(Vec3Caster<Point>(center), _radius, offset.length());
    cylinder->setRotation(osg::Quat(angle, rotation));

    osg::Geode* geode = new osg::Geode;
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(cylinder);
    osg::Vec4 transparentColor = color;
    transparentColor.a() = 0.6f;
    drawable->setColor(transparentColor);
    geode->addDrawable(drawable);

    group->addChild(geode);

    return group;
}

void Cylinder::prepareParameters()
{
    setParameter(0, _normal.x());
    setParameter(1, _normal.y());
    setParameter(2, _normal.z());
    setParameter(3, _point.x());
    setParameter(4, _point.y());
    setParameter(5, _point.z());
    setParameter(6, _radius);
    setParameter(7, 0.0);
}

void Cylinder::applyParameters()
{
    _normal = Vector(getParameter(0), getParameter(1), getParameter(2));
    _point = Point(getParameter(3), getParameter(4), getParameter(5));
    _radius = getParameter(6);
}

bool Cylinder::loadParameters(const std::string& line)
{
    std::stringstream sin(line.substr(line.find_first_of(" \t")+1));
    sin >> _normal >> _point >> _radius;
    return true;
}

bool Cylinder::saveParameters(std::ofstream& fout) const
{
    fout << "# cylinder normal_x normal_y normal_z point_x point_y point_z radius" << std::endl;
    fout << "cylinder " << _normal << " " << _point << " " << _radius << std::endl;
    return fout.good();
}

void Cylinder::computePrecision()
{
    Kernel::Line_3 axis(_point, _normal);

    double max = std::numeric_limits<double>::min();
    double min = std::numeric_limits<double>::max();
    for (size_t i = 0, iEnd = _vecPointIdx.size(); i < iEnd; ++ i) {
        double d = CGAL::squared_distance(axis, _vecPointSet[_vecPointIdx[i]]->point);
        max = max < d? d:max;
        min = min > d? d:min;
    }

    _precision = std::max(std::abs(std::sqrt(max)-_radius), std::abs(std::sqrt(min)-_radius));
}