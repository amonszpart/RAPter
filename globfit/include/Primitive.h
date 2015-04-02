#ifndef Primitive_H
#define Primitive_H

#include <vector>
#include <string>
#include <fstream>

#include <osg/Group>

#include "Types.h"

namespace osg {
    class Group;
}

class Primitive
{
public:
    enum PrimitiveType {
        PT_PLANE,
        PT_SPHERE,
        PT_CYLINDER,
        PT_CONE
    };


    Primitive(const std::vector<RichPoint*>& vecPointSet, PrimitiveType primitiveType);
    ~Primitive(void);

    bool load(std::ifstream& fin);
    bool save(std::ofstream& fout) const;

    const std::vector<size_t>&  getPointIdx() const {return _vecPointIdx;}
    double getPrecision() const {return _precision;}
    PrimitiveType getType() const {return _primitiveType;}

    virtual osg::Group* toGeometry(const osg::Vec4& color) = 0;

    virtual void prepareParameters() = 0;
    virtual void applyParameters() = 0;
    /*
    Parameter format:
    Plane:      normalX normalY normalZ empty   empty   empty   distance
    Cylinder:   normalX normalY normalZ pointX  pointY  pointZ  radius
    Cone:       normalX normalY normalZ apexX   apexY   apexZ   angle
    Sphere:     empty   empty   empty   centerX centerY centerZ radius
    */

    virtual bool getNormal(Vector& normal) const {return false;}
    virtual bool getDistance(double& distance) const {return false;}
    virtual bool getRadius(double& radius) const {return false;}
    virtual bool getCenter(Point& center) const {return false;}

    static size_t getNumParameter() {return 8;}
    double getParameter(size_t idx) const {return _parameters[idx];}
    void setParameter(size_t idx, double value) {_parameters[idx] = value;}

    size_t getIdx() const {return _idx;}
    void setIdx(size_t idx) {_idx = idx;}
protected:
    virtual bool loadParameters(const std::string& line) = 0;
    virtual bool saveParameters(std::ofstream& fout) const = 0;
    virtual void computePrecision() = 0;

    const std::vector<RichPoint*>&  _vecPointSet;
    std::vector<size_t>             _vecPointIdx;
    double                          _precision;
    PrimitiveType                   _primitiveType;
    double                          _parameters[8];
    size_t                          _idx;
private:
};

#endif // Primitive_H
