#include <cassert>
#include "RelationVertex.h"

RelationVertex::RelationVertex(size_t idx, size_t primitiveIdx1)
{
    _relationVertexType = RVT_SINGLE_PRIMITIVE;
    _idx = idx;
    _primitiveIdx1 = primitiveIdx1;
    _parent = idx;
}


RelationVertex::RelationVertex(size_t idx, size_t primitiveIdx1, size_t primitiveIdx2)
{
    _relationVertexType = RVT_DOUBLE_PRIMITIVE;
    _idx = idx;
    _primitiveIdx1 = primitiveIdx1;
    _primitiveIdx2 = primitiveIdx2;
    _parent = idx;
}

RelationVertex::~RelationVertex(void)
{
}

size_t RelationVertex::getPrimitiveIdx1() const
{
    assert(_relationVertexType>=RVT_SINGLE_PRIMITIVE);
    return _primitiveIdx1;
}

size_t RelationVertex::getPrimitiveIdx2() const
{
    assert(_relationVertexType>=RVT_DOUBLE_PRIMITIVE);
    return _primitiveIdx2;
}

bool RelationVertex::operator<(const RelationVertex& other) const
{
    return _primitiveIdx1 < other._primitiveIdx1
        || (_primitiveIdx1 == other._primitiveIdx1 && _primitiveIdx2 < other._primitiveIdx2);
}