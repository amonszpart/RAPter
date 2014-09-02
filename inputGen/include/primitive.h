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

    //! \brief Plane or line type.
    enum TYPE {
        LINE_2D,  //! <\brief Represent a 2D line, z coordinate = 0
        PLANE_3D  //! <\brief Represent a 3D plane
    };

private:
    vec  _coord,  //! <\brief coordinates of one point on the line
         _dir;    //! <\brief direction of the line
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
          _uid(uid),
          _did(did==-1 ? uid : did),
          _type(type) { }

    template <typename MatrixDerived1, typename MatrixDerived2>
    inline LinearPrimitive(LinearPrimitive<_Scalar>::TYPE type,
                           const MatrixDerived1& coord, // coordinates of one point on the line
                           const MatrixDerived2& dir,   // direction of the line
                           int uid = LinearPrimitive<_Scalar>::getUID(),
                           int did = -1) // default behavior, did = uid
        : _coord(coord),
          _dir(dir),
          _uid(uid),
          _did(did==-1 ? uid : did),
          _type(type) { }

    //! \brief Line or plane \see #TYPE
    inline LinearPrimitive<_Scalar>::TYPE  type() const { return _type; }
    inline LinearPrimitive<_Scalar>::TYPE& type() { return _type; }

    inline vec& coord() { return _coord; }
    inline vec& dir()   { return _dir; }
    inline const vec& coord() const { return _coord; }
    inline const vec& dir()   const { return _dir; }

    inline uint uid() const { return _uid; }

    inline int& did()       { return _did; }
    inline int  did() const { return _did; }

    // Test if
    bool areCompatible(const LinearPrimitive<_Scalar>& other){

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
