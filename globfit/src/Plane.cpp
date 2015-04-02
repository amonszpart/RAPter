#include <sstream>

#include <osg/Geode>
#include <osg/Geometry>

#include <CGAL/convex_hull_2.h>

#include "Plane.h"


Plane::Plane(const std::vector<RichPoint*>& vecPointSet) : Primitive(vecPointSet, PT_PLANE)
{
}


Plane::~Plane(void)
{
}

osg::Group* Plane::toGeometry(const osg::Vec4& color)
{
    osg::Group* group = new osg::Group;
    if (_vecPointIdx.size() < 3) {
        return group;
    }

    Kernel::Plane_3 plane(_normal.x(), _normal.y(), _normal.z(), _distance);
    Point origin = plane.projection(_vecPointSet[_vecPointIdx[0]]->point);

    Vector base1 = plane.base1();
    Vector base2 = plane.base2();
    base1 = base1/std::sqrt(base1.squared_length());
    base2 = base2/std::sqrt(base2.squared_length());

    Kernel::Line_3 baseLine1(origin, base1);
    Kernel::Line_3 baseLine2(origin, base2);

    std::vector<Kernel::Point_2> vec2DCoord;
    for (size_t i = 0, iEnd = _vecPointIdx.size(); i < iEnd; ++ i) {
        Point point = plane.projection(_vecPointSet[_vecPointIdx[i]]->point);
        Vector xVector(origin, baseLine1.projection(point));
        Vector yVector(origin, baseLine2.projection(point));
        double x = std::sqrt(xVector.squared_length());
        double y = std::sqrt(yVector.squared_length());
        x = xVector*base1 < 0? -x:x;
        y = yVector*base2 < 0? -y:y;
        vec2DCoord.push_back(Kernel::Point_2(x, y));
    }

    std::vector<Kernel::Point_2> vecConvexHull;
    CGAL::convex_hull_2(vec2DCoord.begin(), vec2DCoord.end(), std::back_inserter(vecConvexHull));

    osg::Vec3Array* vertices = new osg::Vec3Array;
    for (size_t i = 0, iEnd = vecConvexHull.size(); i < iEnd; ++ i) {
        double x = vecConvexHull[i].x();
        double y = vecConvexHull[i].y();
        Point point = origin + x*base1 + y*base2;

        vertices->push_back(Vec3Caster<Point>(point));
    }


    osg::Vec4Array* colors = new osg::Vec4Array;
    osg::Vec4 transparentColor = color;
    transparentColor.a() = 0.6f;
    colors->push_back(transparentColor);

    osg::Geometry* geometry = new osg::Geometry;
    geometry->setUseDisplayList(true);
    geometry->setUseVertexBufferObjects(true);
    geometry->setVertexArray(vertices);
    geometry->setColorArray(colors);
    geometry->setColorBinding(osg::Geometry::BIND_OVERALL);
    geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POLYGON, 0, vertices->size()));

    osg::Geode* geode = new osg::Geode;
    geode->addDrawable(geometry);
    group->addChild(geode);

    return group;
}

void Plane::prepareParameters()
{
    setParameter(0, _normal.x());
    setParameter(1, _normal.y());
    setParameter(2, _normal.z());
    setParameter(3, 0.0);
    setParameter(4, 0.0);
    setParameter(5, 0.0);
    setParameter(6, _distance);
    setParameter(7, 0.0);
}

void Plane::applyParameters()
{
    _normal = Vector(getParameter(0), getParameter(1), getParameter(2));
    _distance = getParameter(6);
}

bool Plane::loadParameters(const std::string& line)
{
    std::stringstream sin(line.substr(line.find_first_of(" \t")+1));
    sin >> _normal >> _distance;
    return true;
}

bool Plane::saveParameters(std::ofstream& fout) const
{
    fout << "# plane normal_x normal_y normal_z d" << std::endl;
    fout << "plane " << _normal << " " << _distance << std::endl;
    return fout.good();
}

void Plane::computePrecision()
{
    Kernel::Plane_3 plane(_normal.x(), _normal.y(), _normal.z(), _distance);
    double max = std::numeric_limits<double>::min();
    for (size_t i = 0, iEnd = _vecPointIdx.size(); i < iEnd; ++ i) {
        double d = CGAL::squared_distance(plane, _vecPointSet[_vecPointIdx[i]]->point);
        max = max<d ? d:max;
    }

    _precision = std::sqrt(max);
}