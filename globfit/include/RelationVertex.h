#ifndef RelationVertex_H
#define RelationVertex_H

#include <cstddef>

class RelationVertex
{
public:
    RelationVertex() {}
    RelationVertex(size_t idx, size_t primitiveIdx1);
    RelationVertex(size_t idx, size_t primitiveIdx1, size_t primitiveIdx2);
    ~RelationVertex(void);

    size_t getPrimitiveIdx1() const;
    size_t getPrimitiveIdx2() const;

    size_t getParent() const {return _parent;}
    void setParent(size_t parent) {_parent = parent;}

    size_t getIdx() const {return _idx;}
    void setIdx(size_t idx) {_idx = idx;}

    bool operator<(const RelationVertex& other) const;
private:
    enum RelationVertexType {
        RVT_SINGLE_PRIMITIVE,
        RVT_DOUBLE_PRIMITIVE
    };
    RelationVertexType _relationVertexType;

    size_t _idx;
    size_t _parent;
    size_t _primitiveIdx1;
    size_t _primitiveIdx2;
};

#endif // RelationVertex_H
