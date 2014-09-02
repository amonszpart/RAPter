#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include "eigen3/Eigen/Core"

namespace InputGen{

/*!
 * \brief Plane or line primitive.
 *
 * Type is set at construction time
 */
template <typename _Scalar>
class LinearPrimitive{

public:
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, 3, 1> vec;
    typedef Eigen::Matrix<Scalar, 2, 1> vec2;

    //! \brief Plane or line type.
    enum TYPE {
        LINE_2D,  //! <\brief Represent a 2D line, z coordinate = 0
        PLANE_3D  //! <\brief Represent a 3D plane
    };

private:
    vec  _coord,  //! <\brief coordinates of one point on the line
         _dir;    //! <\brief direction of the line
    vec2 _dim;    //! <\brief dimensions (width/length) of the primitive

    uint _uid;    //! <\brief unique identifier
    int  _did;    //! <\brief direction id, can be shared with other lines
    LinearPrimitive<_Scalar>::TYPE _type; //! <\brief Current instance type

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline LinearPrimitive(LinearPrimitive<_Scalar>::TYPE type = LinearPrimitive<_Scalar>::LINE_2D,
                           int uid = LinearPrimitive<_Scalar>::getUID(),
                           int did = -1) // default behavior, did = uid
        : _coord(vec::Zero()),
          _dir(vec::Zero()),
          _dim(vec2::Zero()),
          _uid(uid),
          _did(did==-1 ? uid : did),
          _type(type) { }

    //! \brief Line or plane \see #TYPE
    inline LinearPrimitive<_Scalar>::TYPE  type() const { return _type; }
    inline LinearPrimitive<_Scalar>::TYPE& type() { return _type; }

    //inline vec& coord() { return _coord; }
    //inline vec& dir()   { return _dir; }
    inline const vec&  coord() const { return _coord; }
    inline const vec&  dir()   const { return _dir; }
    inline const vec2& dim()   const { return _dim; }

    template <class vec3Derived>
    inline void setCoord(const vec3Derived& coord )
    { _coord = coord; if (_type == LINE_2D) _coord(2) = 0; }

    template <class vec3Derived>
    inline void setDir  (const vec3Derived& dir )
    { _dir = dir; if (_type == LINE_2D) _dir(2) = 0; }

    template <class vec2Derived>
    inline void setDim  (const vec2Derived& dim )
    { _dim = dim; if (_type == LINE_2D) _dir(1) = 0; }

    inline uint uid() const { return _uid; }

    inline int& did()       { return _did; }
    inline int  did() const { return _did; }

    // Test if
    bool areCompatible(const LinearPrimitive<_Scalar>& other){

    }

    template<template <typename> class DisplayFunctor >
    void displayAsLine() const {
        DisplayFunctor<Scalar>::displayVertex(_coord.data());
        //DisplayFunctor<Scalar>::displayVertex((_coord+0.1*_dir).eval().data());
        DisplayFunctor<Scalar>::displayVertex((_coord+vec (_dir(1), -_dir(0), 0)*_dim(0)).eval().data());
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
