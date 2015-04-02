#ifndef Plane_H
#define Plane_H


#include "Primitive.h"
#include "Types.h"


class Plane : public Primitive
{
public:
    Plane(const std::vector<RichPoint*>& vecPointSet);
    ~Plane(void);

    virtual osg::Group* toGeometry(const osg::Vec4& color);

    virtual void prepareParameters();
    virtual void applyParameters();

    virtual bool getNormal(Vector& normal) const {normal = _normal; return true;}
    virtual bool getDistance(double& distance) const {distance = _distance; return true;}
protected:
    virtual bool loadParameters(const std::string& line);
    virtual bool saveParameters(std::ofstream& fout) const;
    virtual void computePrecision();

private:
    Vector      _normal;
    double      _distance;
};

#endif // Plane_H