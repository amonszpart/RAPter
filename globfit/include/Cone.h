 #ifndef Cone_H
#define Cone_H


#include "Primitive.h"
#include "Types.h"

class Cone : public Primitive
{
public:
    Cone(const std::vector<RichPoint*>& vecPointSet);
    ~Cone(void);

    virtual osg::Group* toGeometry(const osg::Vec4& color);

    virtual void prepareParameters();
    virtual void applyParameters();

    virtual bool getNormal(Vector& normal) const {normal = _normal; return true;}
    virtual bool getCenter(Point& center) const {center = _apex; return true;}
protected:
    virtual bool loadParameters(const std::string& line);
    virtual bool saveParameters(std::ofstream& fout) const;
    virtual void computePrecision();

private:
    Vector      _normal;
    Point       _apex;
    double      _angle;
};

#endif // Cone_H