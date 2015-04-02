#ifndef RelationEdge_H
#define RelationEdge_H

#include <ostream>
#include <vector>
#include "RelationVertex.h"

class Primitive;

class RelationEdge
{
    friend std::ostream& operator<<(std::ostream &out, const RelationEdge& relationEdge);
public:
    enum RelationEdgeType {
        // Orientation
        RET_PARALLEL,
        RET_ORTHOGONAL,
        RET_EQUAL_ANGLE,

        // Placement
        RET_COAXIAL,
        RET_COPLANAR,

        // Equality
        RET_EQUAL_LENGTH,
        RET_EQUAL_RADIUS
    };
    RelationEdge() {}
    RelationEdge(RelationEdgeType relationEdgeType, RelationVertex source, RelationVertex target, double score=1.0):
        _relationEdgeType(relationEdgeType), _source(source), _target(target), _score(score) {}
    ~RelationEdge() {}

    static size_t getNumParameter() {return 8;}
    void dumpData(int* p, int nRowNum, size_t row);

    RelationEdgeType getType() const {return _relationEdgeType;}
    const RelationVertex& getSource() const {return _source;}
    const RelationVertex& getTarget() const {return _target;}
    double getScore() const {return _score;}
    void setScore(double score) {_score = score;}

    bool operator<(const RelationEdge& other) const {return _score < other.getScore();}

private:
    RelationEdgeType    _relationEdgeType;
    RelationVertex      _source;
    RelationVertex      _target;
    double              _score;
};

#endif // RelationEdge_H
