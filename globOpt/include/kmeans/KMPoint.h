#ifndef KMPOINT_H
#define KMPOINT_H

#include <vector>

namespace am
{

    // point with primitive id
    template <class _3DPoint>
    struct KMPoint
    {
            typedef std::vector<KMPoint<_3DPoint> > Vec;

            _3DPoint pnt;
            int     primitive_id;
    };

}

#endif // KMPOINT_H
