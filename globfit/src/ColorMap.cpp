#include "ColorMap.h"

static const double jet[256][3] = {
  {  0,         0,         0.5625,  }, 
  {  0,         0,         0.6250,  }, 
  {  0,         0,         0.6875,  }, 
  {  0,         0,         0.7500,  }, 
  {  0,         0,         0.8125,  }, 
  {  0,         0,         0.8750,  }, 
  {  0,         0,         0.9375,  }, 
  {  0,         0,         1.0000,  }, 
  {  0,         0.0625,    1.0000,  }, 
  {  0,         0.1250,    1.0000,  }, 
  {  0,         0.1875,    1.0000,  }, 
  {  0,         0.2500,    1.0000,  }, 
  {  0,         0.3125,    1.0000,  }, 
  {  0,         0.3750,    1.0000,  }, 
  {  0,         0.4375,    1.0000,  }, 
  {  0,         0.5000,    1.0000,  }, 
  {  0,         0.5625,    1.0000,  }, 
  {  0,         0.6250,    1.0000,  }, 
  {  0,         0.6875,    1.0000,  }, 
  {  0,         0.7500,    1.0000,  }, 
  {  0,         0.8125,    1.0000,  }, 
  {  0,         0.8750,    1.0000,  }, 
  {  0,         0.9375,    1.0000,  }, 
  {  0,         1.0000,    1.0000,  }, 
  {  0.0625,    1.0000,    0.9375,  }, 
  {  0.1250,    1.0000,    0.8750,  }, 
  {  0.1875,    1.0000,    0.8125,  }, 
  {  0.2500,    1.0000,    0.7500,  }, 
  {  0.3125,    1.0000,    0.6875,  }, 
  {  0.3750,    1.0000,    0.6250,  }, 
  {  0.4375,    1.0000,    0.5625,  }, 
  {  0.5000,    1.0000,    0.5000,  }, 
  {  0.5625,    1.0000,    0.4375,  }, 
  {  0.6250,    1.0000,    0.3750,  }, 
  {  0.6875,    1.0000,    0.3125,  }, 
  {  0.7500,    1.0000,    0.2500,  }, 
  {  0.8125,    1.0000,    0.1875,  }, 
  {  0.8750,    1.0000,    0.1250,  }, 
  {  0.9375,    1.0000,    0.0625,  }, 
  {  1.0000,    1.0000,    0,  }, 
  {  1.0000,    0.9375,    0,  }, 
  {  1.0000,    0.8750,    0,  }, 
  {  1.0000,    0.8125,    0,  }, 
  {  1.0000,    0.7500,    0,  }, 
  {  1.0000,    0.6875,    0,  }, 
  {  1.0000,    0.6250,    0,  }, 
  {  1.0000,    0.5625,    0,  }, 
  {  1.0000,    0.5000,    0,  }, 
  {  1.0000,    0.4375,    0,  }, 
  {  1.0000,    0.3750,    0,  }, 
  {  1.0000,    0.3125,    0,  }, 
  {  1.0000,    0.2500,    0,  }, 
  {  1.0000,    0.1875,    0,  }, 
  {  1.0000,    0.1250,    0,  }, 
  {  1.0000,    0.0625,    0,  }, 
  {  1.0000,    0,         0,  }, 
  {  0.9375,    0,         0,  }, 
  {  0.8750,    0,         0,  }, 
  {  0.8125,    0,         0,  }, 
  {  0.7500,    0,         0,  }, 
  {  0.6875,    0,         0,  }, 
  {  0.6250,    0,         0,  }, 
  {  0.5625,    0,         0,  }, 
  {  0.5000,    0,         0  }
};

ColorMap::ColorMap()
{
    for (size_t i = 0, iEnd = 64; i < iEnd; ++ i) {
        _jet.push_back(osg::Vec4(jet[i][0], jet[i][1], jet[i][2], 1.0f));
    }

    _lightBlue = osg::Vec4(0.75f, 0.94f, 1.0f, 1.0f);
    return;
}

const osg::Vec4& ColorMap::getColor(NamedColor nameColor)
{
    switch(nameColor) {
    case(LIGHT_BLUE):
        return _lightBlue;
        break;
    default:
        return _lightBlue;
        break;
    }

    return _lightBlue;
}

const osg::Vec4& ColorMap::getColor(Map map, float value, float low, float high)
{
    switch(map) {
    case(JET):
        return getColor(_jet, value, low, high);
        break;
    default:
        return getColor(_jet, value, low, high);
        break;
    }

    return getColor(_jet, value, low, high);
}

const osg::Vec4& ColorMap::getColor(const std::vector<osg::Vec4>& map, float value, float low, float high)
{
    int numColor = map.size();

    if (low > high) {
        std::swap(low, high);
    }

    int index = std::fabs((value-low)/(high-low)*numColor);
    index = index%numColor;

    return map[index];
}
