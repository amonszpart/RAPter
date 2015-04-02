#define _USE_MATH_DEFINES
#include <cmath>

#include <sstream>

#include <osg/Geode>
#include <osg/Geometry>

#include "Cone.h"


Cone::Cone(const std::vector<RichPoint*>& vecPointSet) : Primitive(vecPointSet, PT_CONE)
{
}


Cone::~Cone(void)
{
}

static double degree_to_radian(double deg) {
    return (M_PI * deg / 180.0f); 
}

static Point rotate_in_radian(const Point& p, const Kernel::Line_3& axis, double angle) {
    Point p1 = axis.projection(p);

    // compute the translation
    Kernel::Aff_transformation_3 mat1(CGAL::Translation(), Point(0.0f, 0.0f, 0.0f) - p1);

    // compute the rotation
    Vector u = axis.to_vector();
    u = u/sqrt(u.squared_length());

    double cos_angle = std::cos(angle);
    double sin_angle = std::sin(angle);

    double rot[3][3];
    rot[0][0] = u.x() * u.x() * (1.0f - cos_angle) + cos_angle;
    rot[0][1] = u.x() * u.y() * (1.0f - cos_angle) - u.z() * sin_angle;
    rot[0][2] = u.x() * u.z() * (1.0f - cos_angle) + u.y() * sin_angle;
    rot[1][0] = u.y() * u.x() * (1.0f - cos_angle) + u.z() * sin_angle;
    rot[1][1] = u.y() * u.y() * (1.0f - cos_angle) + cos_angle;
    rot[1][2] = u.y() * u.z() * (1.0f - cos_angle) - u.x() * sin_angle;
    rot[2][0] = u.z() * u.x() * (1.0f - cos_angle) - u.y() * sin_angle;
    rot[2][1] = u.z() * u.y() * (1.0f - cos_angle) + u.x() * sin_angle;
    rot[2][2] = u.z() * u.z() * (1.0f - cos_angle) + cos_angle;

    Kernel::Aff_transformation_3 mat2( 
        rot[0][0], rot[0][1], rot[0][2],  
        rot[1][0], rot[1][1], rot[1][2], 
        rot[2][0], rot[2][1], rot[2][2] );

    Kernel::Aff_transformation_3 mat3( CGAL::Translation(), p1 - Point(0.0f, 0.0f, 0.0f));

    Kernel::Aff_transformation_3 mat = mat3 * mat2 * mat1;
    return mat.transform(p);
}


static Point rotate_in_degree(const Point& p, const Kernel::Line_3& axis, double angle) {
    double rad = degree_to_radian(angle);  
    return rotate_in_radian(p, axis, rad);
}

static osg::Group* makeConialFrustum(Point top, Point base, double topRadius, double baseRadius, osg::Vec4 color)
{
    osg::Group* group = new osg::Group;

    Vector normal(top, base);
    Kernel::Line_3  axis(top, normal);
    Kernel::Plane_3 plane(top, normal);
    Vector perp = plane.base1();
    perp = perp/sqrt(perp.squared_length());
    Point baseSeed = top + perp * topRadius;
    Point topSeed = base + perp * baseRadius;

    osg::Vec3Array* topDisk = new osg::Vec3Array;
    osg::Vec3Array* baseDisk = new osg::Vec3Array;
    osg::Vec3Array* frustum = new osg::Vec3Array;

    unsigned int count = 36; 
    double ang_part = 360.0 / count;
    Point basePrevious = rotate_in_degree(baseSeed, axis, -ang_part);
    Point topPrevious = rotate_in_degree(topSeed, axis, -ang_part);
    for (unsigned int i = 0; i < count; ++ i) {
        Point basePoint = rotate_in_degree(baseSeed, axis, i*ang_part);
        Point topPoint = rotate_in_degree(topSeed, axis, i*ang_part);

        topDisk->push_back(Vec3Caster<Point>(topPoint));
        baseDisk->push_back(Vec3Caster<Point>(basePoint));

        frustum->push_back(Vec3Caster<Point>(topPoint));
        frustum->push_back(Vec3Caster<Point>(basePoint));
        frustum->push_back(Vec3Caster<Point>(basePrevious));
        frustum->push_back(Vec3Caster<Point>(topPrevious));

        basePrevious = basePoint;
        topPrevious = topPoint;
    }

    osg::Vec4Array* colorArray = new osg::Vec4Array;
    colorArray->push_back(color);

    osg::Geode* topGeode = new osg::Geode;
    osg::Geometry* topGeometry = new osg::Geometry;
    topGeometry->setVertexArray(topDisk);
    topGeometry->setColorArray(colorArray);
    topGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
    topGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POLYGON, 0, topDisk->size()));
    topGeode->addDrawable(topGeometry);

    osg::Geode* baseGeode = new osg::Geode;
    osg::Geometry* baseGeometry = new osg::Geometry;
    baseGeometry->setVertexArray(baseDisk);
    baseGeometry->setColorArray(colorArray);
    baseGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
    baseGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POLYGON, 0, baseDisk->size()));
    baseGeode->addDrawable(baseGeometry);

    osg::Geode* frustumGeode = new osg::Geode;
    osg::Geometry* frustumGeometry = new osg::Geometry;
    frustumGeometry->setVertexArray(frustum);
    frustumGeometry->setColorArray(colorArray);
    frustumGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
    frustumGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUADS, 0, frustum->size()));
    frustumGeode->addDrawable(frustumGeometry);

    group->addChild(topGeode);
    group->addChild(baseGeode);
    group->addChild(frustumGeode);

    return group;
}

osg::Group* Cone::toGeometry(const osg::Vec4& color)
{
    osg::Group* group = new osg::Group;
    if (_vecPointIdx.size() < 2) {
        return group;
    }

    Kernel::Line_3 axis(_apex, _normal);
    double topHeight = std::numeric_limits<double>::max();
    double baseHeight = 0;
    for (size_t i = 0, iEnd = _vecPointIdx.size(); i < iEnd; ++ i) {
        Point projection = axis.projection(_vecPointSet[_vecPointIdx[i]]->point);
        Vector vector(_apex, projection);
        double d = vector*_normal;
        if (abs(d) < abs(topHeight)) {
            topHeight = d;
        }
        if (abs(d) > abs(baseHeight)) {
            baseHeight = d;
        }
    }

    Point top = _apex + topHeight*_normal;
    Point base = _apex + baseHeight*_normal;

    osg::Vec4 transparentColor = color;
    transparentColor.a() = 0.6f;

    return makeConialFrustum(top, base, topHeight*tan(_angle), baseHeight*tan(_angle), transparentColor);
}

void Cone::prepareParameters()
{
    setParameter(0, _normal.x());
    setParameter(1, _normal.y());
    setParameter(2, _normal.z());
    setParameter(3, _apex.x());
    setParameter(4, _apex.y());
    setParameter(5, _apex.z());
    setParameter(6, _angle);
    setParameter(7, 0.0);
}

void Cone::applyParameters()
{
    _normal = Vector(getParameter(0), getParameter(1), getParameter(2));
    _apex = Point(getParameter(3), getParameter(4), getParameter(5));
    _angle = getParameter(6);
}

bool Cone::loadParameters(const std::string& line)
{
    std::stringstream sin(line.substr(line.find_first_of(" \t")+1));
    sin >> _normal >> _apex >> _angle;
    return true;
}

bool Cone::saveParameters(std::ofstream& fout) const
{
    fout << "# cone normal_x normal_y normal_z apex_x apex_y apex_z angle" << std::endl;
    fout << "cone " << _normal << " " << _apex << " " << _angle << std::endl;
    return fout.good();
}

void Cone::computePrecision()
{
    Kernel::Line_3 axis(_apex, _normal);
    double cos_theta = cos(_angle/2);
    double sin_theta = sin(_angle/2);

    double max = std::numeric_limits<double>::min();
    for (size_t i = 0, iEnd = _vecPointIdx.size(); i < iEnd; ++ i) {
        const Point& point = _vecPointSet[_vecPointIdx[i]]->point;
        Point projection = axis.projection(point);
        double a = std::sqrt(CGAL::squared_distance(_apex, projection));
        double b = std::sqrt(CGAL::squared_distance(point, projection));
        double d = std::abs(a*sin_theta-b*cos_theta);
        max = max < d? d:max;
    }

    _precision = max;
}