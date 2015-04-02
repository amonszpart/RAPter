#ifndef ColorMap_H
#define ColorMap_H

#include <osg/Vec4>
#include <vector>
#include <cstddef>
#include <cmath>

class ColorMap {
public:
    static ColorMap& Instance() {
        static ColorMap theSingleton;
        return theSingleton;
    }

    typedef enum {JET} Map;
    typedef enum {LIGHT_BLUE} NamedColor;

    /* more (non-static) functions here */
    const osg::Vec4& getColor(NamedColor nameColor);
    const osg::Vec4& getColor(Map map, float value, float low=0.0f, float high=1.0f);

private:
    ColorMap();                            // ctor hidden
    ColorMap(ColorMap const&){}            // copy ctor hidden
    ColorMap& operator=(ColorMap const&){} // assign op. hidden
    ~ColorMap(){}                          // dtor hidden

    std::vector<osg::Vec4> _jet;
    osg::Vec4              _lightBlue;

    const osg::Vec4& getColor(const std::vector<osg::Vec4>& map, float value, float low, float high);
};

#endif // ColorMap_H
