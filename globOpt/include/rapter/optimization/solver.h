#ifndef RAPTER_SOLVER_H__
#define RAPTER_SOLVER_H__

//#include "qcqpcpp/io/io.h"    // read/writeSparseMatrix
#include "Eigen/Sparse"         // Eigen::SparseMatrix (solve)
#include "rapter/typedefs.h"    // rapter::Scalar (solve)

namespace rapter {

/*! \brief Not used.
 */
struct SolverParams
{
        int n_points = 50;
}; // ... struct SolverParams

class Solver
{
    public:
        static const int DO_RETRY = -11; // signal to tell solveCli to run again

        //typedef Eigen::Matrix<Scalar,3,1>                   Vector;
        typedef Eigen::SparseMatrix<Scalar,Eigen::RowMajor> SparseMatrix;

        template < class _PrimitiveContainerT
                 , class _InnerPrimitiveContainerT
                 , class _PrimitiveT
                 >
        static inline int solve      ( int argc, char** argv );

        /*! \brief Globfit planned. \todo: move to datafit.h. */
        template < class _PrimitiveContainerT
                 , class _InnerPrimitiveContainerT
                 , class _PrimitiveT
                 >
        static inline int datafit    ( int argc, char** argv );

        /*! \brief Unused for now. */
        static inline Eigen::Matrix<rapter::Scalar,3,1> checkSolution( std::vector<Scalar>       const& x
                                                                  , SparseMatrix              const& qo
                                                                  , SparseMatrix              const& Qo
                                                                  , SparseMatrix              const& A
                                                                  , Eigen::Matrix<Scalar,3,1> const& weights );
}; //...class Solver

} //... ns rapter

#include "rapter/optimization/impl/solver.hpp"

#endif // RAPTER_SOLVER_H__
