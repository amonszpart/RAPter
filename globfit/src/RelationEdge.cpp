#include "RelationVertex.h"
#include "Primitive.h"
#include "RelationEdge.h"

void RelationEdge::dumpData(int* p, int nRowNum, size_t row)
{
    p[0*nRowNum+row] = _relationEdgeType;

    if (_relationEdgeType == RET_EQUAL_ANGLE || _relationEdgeType == RET_EQUAL_LENGTH) {
        p[1*nRowNum+row] = (int)_target.getPrimitiveIdx1();
        p[2*nRowNum+row] = (int)_target.getPrimitiveIdx2();
        p[3*nRowNum+row] = (int)_source.getPrimitiveIdx1();
        p[4*nRowNum+row] = (int)_source.getPrimitiveIdx2();
    } else {
        p[1*nRowNum+row] = (int)_target.getPrimitiveIdx1();
        p[2*nRowNum+row] = (int)_source.getPrimitiveIdx1();
        p[3*nRowNum+row] = -1;
        p[4*nRowNum+row] = -1;
    }

    return;
}

std::ostream& operator<<(std::ostream &out, const RelationEdge& relationEdge)
{
    if (relationEdge.getType() == RelationEdge::RET_EQUAL_ANGLE || relationEdge.getType() == RelationEdge::RET_EQUAL_LENGTH) {
        out << relationEdge.getType() << " "
            << relationEdge.getTarget().getPrimitiveIdx1() << " " << relationEdge.getTarget().getPrimitiveIdx2() << " "
            << relationEdge.getSource().getPrimitiveIdx1() << " " << relationEdge.getSource().getPrimitiveIdx2();
    } else {
        out << relationEdge.getType() << " "
            << relationEdge.getTarget().getPrimitiveIdx1() << " " << relationEdge.getSource().getPrimitiveIdx1() << " "
            << "-1 -1";
    }
    return out;
}