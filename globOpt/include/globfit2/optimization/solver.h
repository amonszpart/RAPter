#ifndef __GF2_SOLVER_H__
#define __GF2_SOLVER_H__

//////////////
/// Solver
//////////////

#ifdef GF2_USE_GUROBI
#   include "optimization/qp/gurobiOpt.h"
#endif
//#include "globfit2/primitives/pointPrimitive.h"
//#include "globfit2/primitives/linePrimitive2.h" // remove, if typedef is moved
#include "qcqpcpp/io/io.h"                      // read/writeSparseMatrix
#include "globfit2/globOpt_types.h" // _2d::PrimitiveT, etc.

namespace GF2 {

/*! \brief Not used.
 */
struct SolverParams
{
        int n_points = 50;
}; // ... struct SolverParams

class Solver
{
    public:
        //typedef GF2::Scalar                                 Scalar;             // in globOpt_types.h
        typedef Eigen::Matrix<Scalar,3,1>                   Vector;
        //typedef LinePrimitive2                              PrimitiveT;
        //typedef PointPrimitive                              PointPrimitiveT;
        //typedef std::vector<std::vector<PrimitiveT> >       PrimitiveContainerT;
        //typedef std::vector<PointPrimitiveT>                PointContainerT;
        typedef Eigen::SparseMatrix<Scalar,Eigen::RowMajor> SparseMatrix;

        //static inline int show       ( int argc, char** argv );
#if WITH_SAMPLE_INPUT
        static inline int sampleInput( int argc, char** argv );
#endif // WITH_SAMPLE_INPUT
        template < class    _PrimitiveContainerT
                 , class    _PointContainerT
                 , typename _Scalar
                 , class    _PointPrimitiveT
                 , class    _PrimitiveT
                 >
        static inline int generateCli   ( int argc, char** argv );
        //static inline int formulate  ( int argc, char** argv );
        template < class _PrimitiveContainerT
                 , class _InnerPrimitiveContainerT
                 , class _PrimitiveT
                 >
        static inline int solve      ( int argc, char** argv );
        template < class _PrimitiveContainerT
                 , class _InnerPrimitiveContainerT
                 , class _PrimitiveT
                 >
        static inline int datafit    ( int argc, char** argv );

        //static inline int run        ( std::string img_path, Scalar const scale, std::vector<Scalar> const& angles, int argc, char** argv ) __attribute__ ((deprecated));

        static inline Eigen::Matrix<GF2::Scalar,3,1> checkSolution( std::vector<Scalar>       const& x
                                                                  , SparseMatrix              const& qo
                                                                  , SparseMatrix              const& Qo
                                                                  , SparseMatrix              const& A
                                                                  , Eigen::Matrix<Scalar,3,1> const& weights );
}; //...class Solver

} //... ns GF2

#ifndef GF2_INC_SOLVER_HPP
#   define GF2_INC_SOLVER_HPP
#   include "globfit2/optimization/impl/solver.hpp"
#endif

#endif // __GF2_SOLVER_H__
