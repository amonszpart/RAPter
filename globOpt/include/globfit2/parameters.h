#ifndef GF2_PARAMETERS_H
#define GF2_PARAMETERS_H

#include "Eigen/Dense"
#include "globfit2/optimization/problemSetup.h"

namespace GF2 {

//! \brief Collection of parameters the \ref ProblemSetup::formulate needs.
template <typename _Scalar>
struct ProblemSetupParams
{
    //! \brief Scale parameter of input.
    _Scalar                      scale           = 0.05;
    //! \brief Desired angles.
    std::vector<_Scalar>         angles          = {0, M_PI_2, M_PI };
    //! \brief Optimization weights. [ 0: unary, 1: pairwise, 2: complexity ]
    Eigen::Matrix<_Scalar,-1,1>  weights         = (Eigen::Matrix<_Scalar,-1,1>(3,1) << 10000, 10, 1).finished();
    //! \brief Options defined in \ref ProblemSetup::CONSTR_MODE
    ProblemSetup::CONSTR_MODE    constr_mode     = ProblemSetup::PATCH_WISE;
    //! \brief Optoins defined \ref ProblemSetup::DATA_COST_MODE
    ProblemSetup::DATA_COST_MODE data_cost_mode  = ProblemSetup::ASSOC_BASED;
    //! \brief Point-count threshold, that decides if a patch is big or small.
    int                          patch_pop_limit = 5;
};

template <typename _Scalar>
struct MergeParams
{
    //! \brief Scale parameter of input.
    _Scalar                      scale           = 0.05;
    //! \brief Desired angles.
    std::vector<_Scalar>         angles          = {0, M_PI_2, M_PI };
};

}

#endif // GF2_PARAMETERS_H
