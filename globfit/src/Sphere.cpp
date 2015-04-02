#include <sstream>

#include <osg/Geode>
#include <osg/Shape>
#include <osg/ShapeDrawable>

#include "Types.h"

#include "Sphere.h"


Sphere::Sphere(const std::vector<RichPoint*>& vecPointSet) : Primitive(vecPointSet, PT_SPHERE)
{
}


Sphere::~Sphere(void)
{
}

osg::Group* Sphere::toGeometry(const osg::Vec4& color)
{
    osg::Geode* geode = new osg::Geode;
    osg::Sphere* sphere = new osg::Sphere(Vec3Caster<Point>(_center), _radius);
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(sphere);
    osg::Vec4 transparentColor = color;
    transparentColor.a() = 0.6f;
    drawable->setColor(transparentColor);
    geode->addDrawable(drawable);

    osg::Group* group = new osg::Group;
    group->addChild(geode);

    return group;
}

void Sphere::prepareParameters()
{
    setParameter(0, 0.0);
    setParameter(1, 0.0);
    setParameter(2, 0.0);
    setParameter(3, _center.x());
    setParameter(4, _center.y());
    setParameter(5, _center.z());
    setParameter(6, _radius);
    setParameter(7, 0.0);
}

void Sphere::applyParameters()
{
    _center = Point(getParameter(3), getParameter(4), getParameter(5));
    _radius = getParameter(6);
}

bool Sphere::loadParameters(const std::string& line)
{
    std::stringstream sin(line.substr(line.find_first_of(" \t")+1));
    sin >> _center >> _radius;
    return true;
}

bool Sphere::saveParameters(std::ofstream& fout) const
{
    fout << "# sphere center_x center_y center_z radius" << std::endl;
    fout << "sphere " << _center << " " << _radius << std::endl;
    return fout.good();
}

void Sphere::computePrecision()
{
    double distance = 0;
    for (size_t i = 0, iEnd = _vecPointIdx.size(); i < iEnd; ++ i) {
        double d = CGAL::squared_distance(_center, _vecPointSet[_vecPointIdx[i]]->point);
        distance += std::sqrt(d);
    }

    _precision = std::abs(distance/_vecPointIdx.size()-_radius);
}