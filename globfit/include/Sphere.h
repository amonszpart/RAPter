#ifndef Sphere_H
#define Sphere_H


#include "Primitive.h"
#include "Types.h"


class Sphere : public Primitive
{
public:
    Sphere(const std::vector<RichPoint*>& vecPointSet);
    ~Sphere(void);

    virtual osg::Group* toGeometry(const osg::Vec4& color);

    virtual void prepareParameters();
    virtual void applyParameters();

    virtual bool getRadius(double& radius) const {radius = _radius; return true;}
    virtual bool getCenter(Point& center) const {center = _center; return true;}

protected:
    virtual bool loadParameters(const std::string& line);
    virtual bool saveParameters(std::ofstream& fout) const;
    virtual void computePrecision();

private:
    Point       _center;
    double      _radius;
};

#endif // Sphere_H