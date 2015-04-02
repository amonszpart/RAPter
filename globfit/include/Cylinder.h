#ifndef Cylinder_H
#define Cylinder_H

#include "Primitive.h"
#include "Types.h"

class Cylinder : public Primitive
{
public:
    Cylinder(const std::vector<RichPoint*>& vecPointSet);
    ~Cylinder(void);

    virtual osg::Group* toGeometry(const osg::Vec4& color);

    virtual void prepareParameters();
    virtual void applyParameters();

    virtual bool getNormal(Vector& normal) const {normal = _normal; return true;}
    virtual bool getRadius(double& radius) const {radius = _radius; return true;}
    virtual bool getCenter(Point& center) const {center = _point; return true;}

protected:
    virtual bool loadParameters(const std::string& line);
    virtual bool saveParameters(std::ofstream& fout) const;
    virtual void computePrecision();

private:
    Vector      _normal;
    Point       _point;
    double      _radius;
};

#endif // Cylinder_H