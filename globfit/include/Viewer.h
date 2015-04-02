#ifndef Viewer_H
#define Viewer_H

#include <cmath>

namespace osg {
    class Node;
}

class Viewer
{
public:
    Viewer();
    ~Viewer(void);

    void operator()(std::vector<osg::Node*> vecViewData);
};

#endif // Viewer_H
