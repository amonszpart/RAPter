#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include "Eigen/Core"
#include "convexHull2D.h"
#include <iostream>


namespace InputGen{

namespace internal_primitive{
template <typename _Scalar>
struct Types{
    typedef Eigen::Matrix<_Scalar, 3, 1> primitive_vec;
    typedef Eigen::Matrix<_Scalar, 2, 1> primitive_vec2;
};
}


template <typename _Scalar>
struct PiecewiseLinear{
    typedef typename internal_primitive::Types<_Scalar>::primitive_vec vec;
    virtual bool isInside ( const vec& v ) const = 0;
    virtual vec getNormal() const = 0;
    virtual vec getCenter() const = 0;
    virtual _Scalar getHeight() const = 0; // get height along the normal direction

    virtual bool isCompatible(const PiecewiseLinear<_Scalar>& other) const {
        return std::abs(getNormal().dot(other.getNormal())) >= _Scalar(1) - Eigen::NumTraits<_Scalar>::dummy_precision() &&
                getHeight() == other.getHeight();
    }
};

template <typename _Scalar, class Derived>
struct PiecewiseLinearCRTP : public PiecewiseLinear<_Scalar>{
    typedef _Scalar Scalar;
    typedef typename internal_primitive::Types<_Scalar>::primitive_vec vec;

    virtual bool isInside ( const vec& v ) const
    { return (static_cast<const Derived*>(this))->_isInside(v); }
};

template <typename _Scalar, int _NbPoly>
struct PlanarConvexPolygon :
        public PiecewiseLinearCRTP<_Scalar, PlanarConvexPolygon<_Scalar, _NbPoly> >
{
    enum { NbPoly=_NbPoly };
    typedef PiecewiseLinearCRTP<_Scalar, PlanarConvexPolygon<_Scalar, _NbPoly> > Base;
    typedef typename Base::vec vec;
    typedef _Scalar Scalar;

    template<typename VertexContainer>
    inline PlanarConvexPolygon(int* ids = NULL, const VertexContainer& vertices = std::vector<vec>()): Base() {
        static_assert( NbPoly >= 3, "PlanarPolygon requires at least 3 vertices");

        if (ids != NULL && vertices.size() != 0){
            for (int i = 0; i!= NbPoly; ++i)
                _vertices[i] = vertices[ids[i]];
            recomputeHull();
        }
    }



//    virtual void setVertex(vec& v, int index){
//        if (index < NbPoly) {
//            _vertices[index] = v;
//            recomputeHull();
//        }
//    }

    virtual vec getNormal() const { return _vertices[0].cross(_vertices[1]).normalized(); }
    virtual Scalar getHeight() const { return _gcenter.dot(getNormal()); }
    virtual vec getCenter() const { return _gcenter; }



    typedef ConvexHull2D<typename internal_primitive::Types<_Scalar>::primitive_vec2> ConvexHull;

    bool inline _isInside(const vec &v) const { return _hull.isInside(moveToLocalAndProject(v, getNormal())); }

protected:

    inline
    typename ConvexHull::Point
    moveToLocalAndProject(const vec &p, const vec& n) const {
//        vec proj = (p - _gcenter)       // get q = p expressed in local space
//                .dot(n)              // get orthogonal distance to the plane
//                * -n                 // get dq = proj(q) - q
//                + (p - _gcenter)     // get proj(q)
//                ;
//        std::cout << p.transpose() << " -> " << proj.transpose()
//                  << (_hull.isInside(proj.head(2)) ? " GOOD ! ": "")
//                  << std::endl;
        return ((p - _gcenter)       // get q = p expressed in local space
                .dot(n)              // get orthogonal distance to the plane
                * -n                 // get dq = proj(q) - q
                + (p - _gcenter)     // get proj(q)
                ).head(2).eval();    // extract two first coordinates, the third one being = 0
    }

    void inline recomputeHull(){
        // we check we have different points in the cloud
        if (_vertices[0] != _vertices[1] &&
            _vertices[0] != _vertices[2] ){

            // estimate plane direction using the three first coordinates
            vec n = getNormal();

            // compute centroid
            vec _gcenter (vec::Zero());
            std::for_each(_vertices.begin(), _vertices.end(), [&_gcenter] (const vec& p){ _gcenter+=p; });
            _gcenter /= _Scalar(_vertices.size());

            // iterate over all points, move them in the local frame, project, and add to the convex hull
            std::vector< typename ConvexHull::Point > projectedSet;
            projectedSet.reserve(_vertices.size());
            std::for_each(_vertices.begin(), _vertices.end(), [&_gcenter, &projectedSet, &n, this] (const vec& p)
            {
                projectedSet.push_back( this->moveToLocalAndProject (p, n) );
            });

            _hull.compute(projectedSet);
        }else {
            _gcenter = vec::Zero();
        }
    }

    //! \brief Array containing points, all expected to lie on a plane
    std::array<vec, NbPoly> _vertices;
    vec _gcenter;

    ConvexHull _hull;
};


/*!
 * \brief Plane or line primitive.
 *
 * Type is set at construction time
 */
template <typename _Scalar>
class LinearPrimitive{

public:
    typedef _Scalar Scalar;
    typedef typename internal_primitive::Types<_Scalar>::primitive_vec  vec;
    typedef typename internal_primitive::Types<_Scalar>::primitive_vec2 vec2;

    //! \brief Plane or line type.
    enum TYPE {
        LINE_2D,  //! <\brief Represent a 2D line, z coordinate = 0
        PLANE_3D  //! <\brief Represent a 3D plane
    };

private:
    vec  _coord,  //! <\brief coordinates of one point on the line
         _normal; //! <\brief direction of the line
    vec2 _dim;    //! <\brief dimensions (width/length) of the primitive

    uint _uid;    //! <\brief unique identifier
    int  _did;    //! <\brief direction id, can be shared with other lines
    LinearPrimitive<_Scalar>::TYPE _type; //! <\brief Current instance type

    std::vector<PiecewiseLinear<_Scalar>*> _pieceWiseElements;

    inline void _synchInternalStateFromPWElements(){
        if (_pieceWiseElements.empty()) return;

        // midpoint center
        vec midPoint = vec::Zero();
        std::for_each(_pieceWiseElements.begin(), _pieceWiseElements.end(), [&midPoint] (const PiecewiseLinear<_Scalar>* pwl){
            midPoint+=pwl->getCenter(); // this is not exact, it would be better to take the bbox and then the center
        });
        midPoint /= _Scalar(_pieceWiseElements.size());

        // normal
        _normal = _pieceWiseElements.front()->getNormal();

        // just check: we assume that the scene is in the bounding box, so we just create a dummy size of 1,1.
        _coord = - Scalar(0.5) * getTangentVector() + midPoint;

        _dim << 1,1;

        // local dimensions (using hull)
    }

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline LinearPrimitive(LinearPrimitive<_Scalar>::TYPE type = LinearPrimitive<_Scalar>::LINE_2D,
                           int uid = LinearPrimitive<_Scalar>::getUID(),
                           int did = -1) // default behavior, did = uid
        : _coord(vec::Zero()),
          _normal(vec(1,0,0)),
          _dim(vec2::Zero()),
          _uid(uid),
          _did(did==-1 ? uid : did),
          _type(type) { }

    inline LinearPrimitive(PiecewiseLinear<_Scalar>* element,
                           int uid = LinearPrimitive<_Scalar>::getUID(),
                           int did = -1)
        : _coord(vec::Zero()),
          _normal(vec(1,0,0)),
          _dim(vec2::Zero()),
          _uid(uid),
          _did(did==-1 ? uid : did),
          _type(LinearPrimitive<_Scalar>::PLANE_3D)
    {
        addPiecewiseElement(element);
    }

    //! \brief Line or plane \see #TYPE
    inline LinearPrimitive<_Scalar>::TYPE  type() const { return _type; }
    inline LinearPrimitive<_Scalar>::TYPE& type() { return _type; }

    //inline vec& coord() { return _coord; }
    //inline vec& dir()   { return _dir; }
    inline const vec&  coord()   const { return _coord; }
    inline const vec&  normal()  const { return _normal; }
    inline const vec2& dim()     const { return _dim; }

    template <class vec3Derived>
    inline void setCoord(const vec3Derived& coord )
    { _coord = coord; if (_type == LINE_2D) _coord(2) = 0; }

    template <class vec3Derived>
    inline void setNormal  (const vec3Derived& normal )
    { _normal= normal; if (_type == LINE_2D) _normal(2) = 0; }

    template <class vec2Derived>
    inline void setDim  (const vec2Derived& dim )
    { _dim = dim; if (_dim(1) != 0.) _type = PLANE_3D; }

    inline const uint& uid() const { return _uid; }

    inline int& did()       { return _did; }
    inline int  did() const { return _did; }

    template<template <typename> class DisplayFunctor >
    void displayAsLine() const {
        DisplayFunctor<Scalar>::displayVertex(_coord.data());
        //DisplayFunctor<Scalar>::displayVertex((_coord+0.1*_dir).eval().data());
        DisplayFunctor<Scalar>::displayVertex((_coord+getTangentVector()*_dim(0)).eval().data());
    }

    // check if the current sample v is inside the primitive
    inline bool isInside(const vec &v) const{
        if (_pieceWiseElements.empty()) return true; // we don't care and assume the sampling is ok

        return std::find_if( _pieceWiseElements.cbegin(),
                             _pieceWiseElements.cend(),
                             [ &v ] (const PiecewiseLinear<_Scalar>*e) {
            return e->isInside(v);
        }) != _pieceWiseElements.cend();
    }

    const std::vector<PiecewiseLinear<_Scalar>*>& piecewiseElements() const { return _pieceWiseElements; }
    bool addPiecewiseElement(PiecewiseLinear<_Scalar>* e) {
        bool valid = true;

        if (! _pieceWiseElements.empty()){
            // test against existing primitives, valid only if we have at least one primitive compatible
            valid =  find_if( _pieceWiseElements.cbegin(), _pieceWiseElements.cend(),
                              [&e](const PiecewiseLinear<_Scalar>* ref ) { return ref->isCompatible(*e);})
                    != _pieceWiseElements.end();
        }

        if(valid){
            _pieceWiseElements.push_back(e);
            _synchInternalStateFromPWElements();
        }

        return valid;
    }


    //! \brief Compute tangent vector
    inline vec getTangentVector() const {
        // rotation of the normal vector around z axis
        //if (_type == LINE_2D)
        if(_normal(0) != _normal(1))
            return vec(_normal(1), -_normal(0), 0);
        else
            return vec(_normal(2), 0, -_normal(0));

        //else{
        //    std::cerr << "Invalid request in " << __FILE__ << ":" << __LINE__ << std::endl;
        //    return vec::Zero();
        //}

    }

    inline vec getMidPoint() const {
        return _coord+0.5*getTangentVector()*_dim(0);
    }

    inline vec getEndPoint() const {
        return _coord+getTangentVector()*_dim(0);
    }


    static inline int getUID() {
        static int uid = 0;
        int id =  uid;
        uid++;
        return id;
    }



}; // class Primitive
} // namespace InputGen

#endif // PRIMITIVE_H
