#ifndef BONMINOPT_H
#define BONMINOPT_H

#include "coin/BonTMINLP.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "qcqpcpp/optProblem.h"
#include "coin/BonBonminSetup.hpp" // Bonmin::Algorithm enum


#include <chrono>

namespace qcqpcpp
{

class BonminOptException : public std::runtime_error
{
        using std::runtime_error::runtime_error;
};

template <typename _Scalar>
class BonminOpt : public qcqpcpp::OptProblem<_Scalar>
{
        //typedef Eigen::Matrix<_Scalar,-1,1> VectorX;
        typedef typename qcqpcpp::OptProblem<_Scalar>::VectorX VectorX;
        enum { NeedsToAlign = (sizeof(VectorX)%16)==0 };
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

        typedef typename qcqpcpp::OptProblem<_Scalar>               ParentType;
        typedef typename ParentType::SparseMatrix                   SparseMatrix;

        //! \brief BonminOpt    Default constructor.
        BonminOpt()
            : _algCode  ( Bonmin::Algorithm::B_BB )
            , _nodeLimit( 100 )
            , _maxSolutions( 0 )
            , _printSol ( false )
            , _debug    ( false )
        {}

        //! \brief ~BonminOpt   virtual destructor.
        virtual ~BonminOpt() {}

        /** \name Functions from OptProblem. */
        //@{
        //! \brief update               main Sets up problem tailored to solver. Bonminopt needs jacobian and hessian computation.
        //! \param verbose              Controls logging to console.
        virtual int update( bool verbose = false  );

        //! \brief optimize             Runs optimization. Please call update before.
        //! \param x_out                Output values.
        //! \param objective_sense      Minimization or Maximization.
        virtual int optimize( std::vector<_Scalar> *x_out /* = NULL */, typename ParentType::OBJ_SENSE objecitve_sense /* = OBJ_SENSE::MINIMIZE */ );

        virtual _Scalar getINF   () const override { return std::numeric_limits<_Scalar>::max(); } // DBL_MAX
        virtual int     getOkCode() const override { return 0; }
        //@}

        static inline Bonmin::TMINLP::VariableType getVarTypeCustom( typename ParentType::VAR_TYPE var_type );
        //inline SparseMatrix const   getJacobian            ( VectorX const& x ) const;// { return _jacobian; }
        //inline SparseMatrix const&  getHessian             ()        const;// { return _hessian; }
        inline VectorX      const&  getCachedqo            ()        const { return _qo; }
        inline SparseMatrix const&  getCachedQo            ()        const { return _Qo; }
        inline SparseMatrix const&  getCachedA             ()        const { return _A; }
        inline SparseMatrix const&  getCachedQ             ( int const j ) const { return _Qs.at(j); }
        inline size_t               getCachedQSize         ()        const { return _Qs.size(); }
        inline VectorX      const&  getGradF               ()        const { return _grad_f; };
        inline VectorX           &  getGradF               ()              { return _grad_f; };

        inline void printSolutionAtEndOfAlgorithm          ()              { _printSol = true; }
        inline bool isDebug                                ()        const { return _debug; }
        inline void setDebug                               ( bool const debug ) { _debug = debug; }
        inline bool isPrintSol                             ()        const { return _printSol; }
        inline void setAlgorithm                           ( Bonmin::Algorithm alg ) { _algCode = alg; }
        inline void setNodeLimit                           ( int nodeLimit )         { _nodeLimit = nodeLimit; }
        inline void setMaxSolutions                        ( int maxSolutions )         { _maxSolutions = maxSolutions; }


    protected:
        //SparseMatrix                 _jacobian; //!< \brief Cached Jacobian of linear constraints.
        //SparseMatrix                 _hessian;  //!< \brief Cached Objective+quadConstriants Hessian, lower triangle only!

        VectorX                      _qo;       //!< \brief Cached full linear objective vector.
        SparseMatrix                 _Qo;       //!< \brief Cached full quadratic objective matrix.
        SparseMatrix                 _A;        //!< \brief Cached full linear constraint matrix.
        std::vector<SparseMatrix>    _Qs;       //!< \brief Cached full quadratic constraint matrices.

        Bonmin::Algorithm           _algCode;   //!< \brief Stores the chosen algorihtm code. 0 = B_Bb default.
        int                         _nodeLimit; //!< \brief How many nodes bonmin can explore.
        int                         _maxSolutions;

        VectorX                     _grad_f; //!< \brief Caches eval_grad_f output.
    private:
        bool                        _printSol;  //!< \brief Flag, to print x in the end.
        bool                        _debug;


}; //...class BonminOpt

template <typename _Scalar>
class BonminTMINLP : public Bonmin::TMINLP
{
    public:
        typedef typename qcqpcpp::OptProblem<_Scalar>                   ParentType;
        typedef typename ParentType::VectorX                            VectorX;
        typedef typename ParentType::SparseMatrix                       SparseMatrix;
        typedef          Eigen::Map<const Eigen::Matrix<Ipopt::Number,-1,1> > MatrixMapT;
        typedef          Eigen::Map<      Eigen::Matrix<Ipopt::Number,-1,1> > MatrixNonConstMapT;

        BonminTMINLP( BonminOpt<_Scalar> &delegate )
            : _delegate( delegate )
            , _ones( VectorX(delegate.getVarCount(),1) ) // row_vector
        {
            _ones.setConstant( _Scalar(1.) );
        }

        //! \brief ~BonminOpt   virtual destructor.
        virtual ~BonminTMINLP() {}

        /** \name Overloaded functions specific to a TMINLP.*/
        //@{
        //! \brief get_variables_types  Pass the type of the variables (INTEGER, BINARY, CONTINUOUS) to the optimizer.
        //! \param n                    size of var_types (has to be equal to the number of variables in the problem)
        //! \param var_types            types of the variables (has to be filled by function).
        virtual bool get_variables_types( Ipopt::Index n, VariableType* var_types );

        //! \brief get_variables_linearity      Pass info about linear and nonlinear variables.
        virtual bool get_variables_linearity( Ipopt::Index n, Ipopt::TNLP::LinearityType* var_types);

        //! \brief get_constraints_linearity    Pass the type of the constraints (LINEAR, NON_LINEAR) to the optimizer.
        //! \param m                    Size of const_types (has to be equal to the number of constraints in the problem)
        //! \param const_types          Types of the constraints (has to be filled by function).
        virtual bool get_constraints_linearity( Ipopt::Index m, Ipopt::TNLP::LinearityType* const_types);
        //@}

        /** \name Overloaded functions defining a TNLP.
         * This group of function implement the various elements needed to define and solve a TNLP.
         * They are the same as those in a standard Ipopt NLP problem*/
        //@{
        //! \brief get_nlp_info         Method to pass the main dimensions of the problem to Ipopt.
        //! \param n                    Number of variables in problem.
        //! \param m                    Number of constraints.
        //! \param nnz_jac_g            Number of nonzeroes in Jacobian of constraints system.
        //! \param nnz_h_lag            Number of nonzeroes in Hessian of the Lagrangean.
        //! \param index_style          Indicate wether arrays are numbered from 0 (C-style) or from 1 (Fortran).
        //! \return true                in case of success.
        virtual bool get_nlp_info( Ipopt::Index& n, Ipopt::Index&m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style );

        //! \brief get_bounds_info      Method to pass the bounds on variables and constraints to Ipopt.
        //! \param n                    Size of x_l and x_u (has to be equal to the number of variables in the problem)
        //! \param x_l                  Lower bounds on variables (function should fill it).
        //! \param x_u                  Upper bounds on the variables (function should fill it).
        //! \param m                    Size of g_l and g_u (has to be equal to the number of constraints in the problem).
        //! \param g_l                  Lower bounds of the constraints (function should fill it).
        //! \param g_u                  Upper bounds of the constraints (function should fill it).
        //! \return true                in case of success.
        virtual bool get_bounds_info( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u );

        //! \brief get_starting_point   Method to to pass the starting point for optimization to Ipopt.
        //! \param init_x               Do we initialize primals?
        //! \param x                    Pass starting primal points (function should fill it if init_x is 1).
        //! \param m                    Size of lambda (has to be equal to the number of constraints in the problem).
        //! \param init_lambda          Do we initialize duals of constraints?
        //! \param lambda               Lower bounds of the constraints (function should fill it).
        //! \return true                in case of success.
        virtual bool get_starting_point( Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                         bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                         Ipopt::Index m, bool init_lambda,
                                         Ipopt::Number* lambda );

        //! \brief eval_f               Method which compute the value of the objective function at point x.
        //! \param n                    size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param obj_value            Value of objective in x (has to be computed by the function).
        //! \return true                in case of success.
        virtual bool eval_f( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value );

        //! \brief eval_grad_f          Method which compute the gradient of the objective at a point x.
        //! \param n                    Size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param grad_f               Gradient of objective taken in x (function has to fill it).
        //! \return true                in case of success.
        virtual bool eval_grad_f( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);

        //! \brief eval_g               Method which compute the value of the functions defining the constraints at a point x.
        //! \param n                    Size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param m                    Size of array g (has to be equal to the number of constraints in the problem)
        //! \param grad_f               Values of the constraints (function has to fill it).
        //! \return true                In case of success.
        virtual bool eval_g( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g);

        //! \brief eval_jac_g           Method to compute the Jacobian of the functions defining the constraints.
        //! \brief                      If the parameter values==NULL fill the arrays iCol and jRow which store the position of
        //! \brief                      the non-zero element of the Jacobian.
        //! \brief                      If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
        //! \param n                    Size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param m                    Size of array g (has to be equal to the number of constraints in the problem)
        //! \param grad_f               Values of the constraints (function has to fill it).
        //! \return true                in case of success.
        virtual bool eval_jac_g( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                                 Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
                                 Ipopt::Number* values );

        //! \brief eval_h               Method to compute the Jacobian of the functions defining the constraints.
        //! \brief                      If the parameter values==NULL fill the arrays iCol and jRow which store the position of
        //! \brief                      the non-zero element of the Jacobian.
        //! \brief                      If the paramenter values!=NULL fill values with the non-zero elements of the Jacobian.
        //! \param n                    Size of array x (has to be the number of variables in the problem).
        //! \param x                    Point where to evaluate.
        //! \param new_x                Is this the first time we evaluate functions at this point? (in the present context we don't care).
        //! \param m                    Size of array g (has to be equal to the number of constraints in the problem)
        //! \param grad_f               Values of the constraints (function has to fill it).
        //! \return true                in case of success.
        virtual bool eval_h( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                             Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                             bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                             Ipopt::Index* jCol, Ipopt::Number* values );


        //! \brief finalize_solution    Method called by Ipopt at the end of optimization.
        virtual void finalize_solution( TMINLP::SolverReturn status,
                                        Ipopt::Index n, const Ipopt::Number* x, Ipopt::Number obj_value );
        //@}

        virtual const SosInfo      * sosConstraints() const { return NULL; }
        virtual const BranchingInfo* branchingInfo () const { return NULL; }

    protected:
        BonminOpt<_Scalar>  &_delegate;
        VectorX              _ones;
}; // ...class BonminTMINLP

} //...ns qcqpcpp

//_____________________________________________________________________________________________________________________
// HPP
#include <random>
#include "coin/BonBonminSetup.hpp"
#include "coin/BonCbc.hpp"          // Bab
#include "coin/BonBabSetupBase.hpp" // Bonmin::IntParameter

#define MYDEBUG 1
namespace qcqpcpp
{

template <typename _Scalar> int
BonminOpt<_Scalar>::update( bool verbose /* = false */ )
{
//    _jacobian = this->estimateJacobianOfConstraints();
    VectorX ones( this->getVarCount(), 1 ); ones.setConstant( _Scalar(1.) );
    this->getJacobian( ones );
    if ( isDebug() )
    {
        std::cout<<"[" << __func__ << "]: " << "jacobian ok" << this->_jacobian << std::endl; fflush(stdout);
    }

    this->precalcHessianCoeffs();
    if ( isDebug() )
    {
        std::cout<<"[" << __func__ << "]: " << "hessian ok";
        for ( size_t i = 0; i != this->_hessians.size(); ++i )
            std::cout << "_hessians[" << i << "]:\n" << this->_hessians[i] << std::endl;
        fflush(stdout);
    }

    // NOTE: _qo is a columnvector for convenient multiplication reasons, whilst sparsematrices are rowVectors by convention
    _qo = this->getLinObjectivesMatrix();
    _Qo = this->getQuadraticObjectivesMatrix();
    for ( int row = 0; row != _Qo.outerSize(); ++row )
    {
        for ( typename SparseMatrix::InnerIterator it(_Qo,row); it; ++it )
        {
            if ( it.row() < it.col() )
            {
                std::cerr << "Quadratic objective matrix not lower-triangular!" << std::endl; fflush(stderr);
                throw new BonminOptException("Quadratic objective matrix not lower-triangular!");
            }
        }
    }

    _A  = this->getLinConstraintsMatrix();
    int max_Q_with_Nonzero = 0;
    std::vector<SparseMatrix> Qs;
    for ( size_t j = 0; j != this->getConstraintCount(); ++j )
    {
        Qs.push_back( this->getQuadraticConstraintsMatrix(j) );
        if ( Qs.back().nonZeros() > 0 )
        {
            max_Q_with_Nonzero = j;
            if ( verbose ) std::cout << "nonzeros: " << Qs.back().nonZeros() << std::endl;
        }
    }

    for ( int j = 0; j != max_Q_with_Nonzero+1; ++j )
    {
        _Qs.push_back( Qs[j] );
    }

    this->_updated = true;

    return EXIT_SUCCESS;
}

template <typename _Scalar> int
BonminOpt<_Scalar>::optimize( std::vector<_Scalar> *x_out /* = NULL */, typename ParentType::OBJ_SENSE objective_sense /* = MINIMIZE */ )
{
    if ( !this->_updated )
    {
        std::cerr << "[" << __func__ << "]: " << "Please call update() first!" << std::endl;
        return EXIT_FAILURE;
    }

    // Now initialize from tminlp
    Bonmin::BonminSetup bonmin2;
    bonmin2.initializeOptionsAndJournalist();

    switch ( this->_algCode )
    {
        case Bonmin::Algorithm::B_BB:
            bonmin2.readOptionsString("bonmin.algorithm B-BB\n"); break;
        case Bonmin::Algorithm::B_OA:
            bonmin2.readOptionsString("bonmin.algorithm B-OA\n"); break;
        case Bonmin::Algorithm::B_QG:
            bonmin2.readOptionsString("bonmin.algorithm B-QG\n"); break;
        case Bonmin::Algorithm::B_Hyb:
            bonmin2.readOptionsString("bonmin.algorithm B-Hyb\n"); break;
        case Bonmin::Algorithm::B_Ecp:
            bonmin2.readOptionsString("bonmin.algorithm B-Ecp\n"); break;
        case Bonmin::Algorithm::B_IFP:
            bonmin2.readOptionsString("bonmin.algorithm B-IFP\n"); break;
        default:
            std::cerr << "[" << __func__ << "]: " << "UNCRECOGNIZED algorithm" << std::endl;
            return EXIT_FAILURE;
            break;
    }

    // need this relay, otherwise, we'll end up with a double free corruption thing
    Ipopt::SmartPtr<BonminTMINLP<_Scalar> > problem = new BonminTMINLP<_Scalar>( *this );
    //bonmin2.options()->SetStringValue("derivative_test","second-order");
    //bonmin2.options()->SetStringValue("derivative_test_print_all","yes");

    std::cout << "[" << __func__ << "]: " << "bonmin.initialize( problem ) called..." << std::endl; fflush(stdout);
    bonmin2.initialize( GetRawPtr(problem) );
    std::cout << "[" << __func__ << "]: " << "bonmin.initialize( problem ) finished..." << std::endl; fflush(stdout);
    //bonmin2.setDoubleParameter( Bonmin::BabSetupBase::AllowableFractionGap, 1e-20 );

    if ( this->getTimeLimit() > _Scalar(0) )
        bonmin2.setDoubleParameter( Bonmin::BabSetupBase::MaxTime, this->getTimeLimit() );
    if ( this->getTolRelGap() > _Scalar(0) )
        bonmin2.setDoubleParameter( Bonmin::BabSetupBase::AllowableFractionGap, this->getTolRelGap());
    std::cout << "[" << __func__ << "]: " << "Bonmin::BabLogLevel= " << bonmin2.getIntParameter( Bonmin::BabSetupBase::BabLogLevel) << std::endl;
    bonmin2.setIntParameter( Bonmin::BabSetupBase::BabLogLevel, 5 );
    //bonmin2.setIntParameter( Bonmin::BabSetupBase::MaxNodes, 1000 );
    bonmin2.setIntParameter( Bonmin::BabSetupBase::BabLogInterval, 10 );
    bonmin2.setIntParameter( Bonmin::BabSetupBase::MaxNodes, _nodeLimit );
    if ( _maxSolutions )
        bonmin2.setIntParameter( Bonmin::BabSetupBase::MaxSolutions, _maxSolutions );

    std::cout << "[" << __func__ << "]: " << "Bonmin::MaxNode = " << bonmin2.getIntParameter( Bonmin::BabSetupBase::MaxNodes ) << ", _nodeLimit: " << _nodeLimit << std::endl;
    std::cout << "[" << __func__ << "]: " << "Bonmin::MaxIterations = " << bonmin2.getIntParameter( Bonmin::BabSetupBase::MaxIterations ) << std::endl;
    std::cout << "[" << __func__ << "]: " << "Bonmin::BabLogInterval = " << bonmin2.getIntParameter( Bonmin::BabSetupBase::BabLogInterval ) << std::endl;

    std::cout << "[" << __func__ << "]: " << "Bonmin::CutoffDecr = " << bonmin2.getDoubleParameter( Bonmin::BabSetupBase::CutoffDecr ) << std::endl;
    std::cout << "[" << __func__ << "]: " << "Bonmin::Cutoff = " << bonmin2.getDoubleParameter( Bonmin::BabSetupBase::Cutoff ) << std::endl;
    std::cout << "[" << __func__ << "]: " << "Bonmin::AllowableGap = " << bonmin2.getDoubleParameter( Bonmin::BabSetupBase::AllowableGap ) << std::endl;
    std::cout << "[" << __func__ << "]: " << "Bonmin::AllowableFractionGap = " << bonmin2.getDoubleParameter( Bonmin::BabSetupBase::AllowableFractionGap ) << std::endl;
    std::cout << "[" << __func__ << "]: " << "Bonmin::IntTol = " << bonmin2.getDoubleParameter( Bonmin::BabSetupBase::IntTol ) << std::endl;
    std::cout << "[" << __func__ << "]: " << "Bonmin::MaxTime = " << bonmin2.getDoubleParameter( Bonmin::BabSetupBase::MaxTime ) << std::endl;

    /*  bb_log_interval: [ Interval at which node level output is printed (number of nodes) {100} ]
                       bb_log_level: [ Level of branch and bound log detail (0-5) {1} ]
                       lp_log_level: [ Level of LP solver log detail (0-4) {0} ]
                     milp_log_level: [ Level of MILP solver log detail (0-4) {0} ]
                      nlp_log_level: [ Level of NLP solver log detail (0-2) {0} ]
                    nlp_log_at_root: [ Level of NLP solver log detail at root node (0-12) {0} ]
                   oa_log_frequency: [ Frequency (in seconds) of OA log messages  {100} ]
                       oa_log_level: [ Level of OA decomposition log detail (0-2) {1} ]*/



    //derivative_test_print_all
    //bonmin.readOptionsString("bonmin.algorithm B-BB\n");
    //bonmin.options()->SetNumericValue("bonmin.time_limit", 5); //changes bonmin's time limit
    //bonmin.options()->SetStringValue("mu_oracle","loqo");

    // Set up done, now let's branch and bound
    try
    {
        Bonmin::Bab bb;
        //bb.setUsingCouenne( true ); // testing

        bb( bonmin2 ); // process parameter file using Ipopt and do branch and bound using Cbc

        std::cout << "[" << __func__ << "]: " << "bonmin finished" << std::endl; fflush(stdout);
    }
    catch ( Bonmin::TNLPSolver::UnsolvedError *E )
    {
        // There has been a failure to solve a problem with Ipopt.
        std::cerr << "Ipopt has failed to solve a problem" << std::endl;
    }
    catch ( Bonmin::OsiTMINLPInterface::SimpleError &E )
    {
        std::cerr << E.className() << "::" << E.methodName()
                  << std::endl
                  << E.message()
                  << std::endl;
    }
    catch ( CoinError &E )
    {
        std::cerr << E.className() << "::" << E.methodName()
                  << std::endl
                  << E.message()
                  << std::endl;
    }
    catch ( std::exception &ex )
    {
        std::cerr << "[" << __func__ << "]: " << "exception: " << ex.what() << std::endl;
    }

    // output answer
    if ( x_out )
    {
        x_out->clear();
        x_out->reserve( this->_x.size() );
        for ( int i = 0; i != this->_x.size(); ++i )
            x_out->push_back( this->_x[i] );
    }

    return !( this->_x.size() == static_cast<typename VectorX::Index>(this->getVarCount()) );
}

template <typename _Scalar> Bonmin::TMINLP::VariableType
BonminOpt<_Scalar>::getVarTypeCustom( typename ParentType::VAR_TYPE var_type )
{
    switch ( var_type )
    {
        case ParentType::VAR_TYPE::CONTINUOUS:
            return Bonmin::TMINLP::VariableType::CONTINUOUS;
            break;
        case ParentType::VAR_TYPE::INTEGER:
            return Bonmin::TMINLP::VariableType::INTEGER;
            break;
        case ParentType::VAR_TYPE::BINARY:
            return Bonmin::TMINLP::VariableType::BINARY;
            break;
        default:
            std::cerr << "[" << __func__ << "]: " << "Unrecognized file type, returning continuous" << std::endl;
            return Bonmin::TMINLP::VariableType::CONTINUOUS;
            break;
    } //... switch

} //...BonminOpt::getVarTypeCustom()

//_____________________________________________________________________________________________________________________

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_variables_types( Ipopt::Index n, VariableType* var_types )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n << ")" << std::endl; fflush( stdout );
    }

    if ( n != static_cast<Ipopt::Index>(_delegate.getVarCount()) )
        throw new BonminOptException( "[BonminOpt::get_variables_types] n != getVarCount()" );

    for ( size_t j = 0; j != _delegate.getVarCount(); ++j )
        var_types[ j ] = BonminOpt<_Scalar>::getVarTypeCustom( _delegate.getVarType(j) );

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_variables_type()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_variables_linearity( Ipopt::Index n, Ipopt::TNLP::LinearityType* var_types )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n << ")" << std::endl; fflush( stdout );
    }

    if ( n != static_cast<Ipopt::Index>(_delegate.getVarCount()) )
        throw new BonminOptException( "[BonminOpt::get_variables_types] n != getVarCount()" );

    // NOTE: this is hard coded for now, non-lin cases not tested
    for ( size_t j = 0; j != _delegate.getVarCount(); ++j )
        var_types[j] = Ipopt::TNLP::NON_LINEAR;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_variables_linearity()


template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_constraints_linearity( Ipopt::Index m, Ipopt::TNLP::LinearityType* const_types )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(m = " << m << ")" << std::endl; fflush( stdout );
    }

    if ( m != static_cast<Ipopt::Index>(_delegate.getConstraintCount()) )
        throw new BonminOptException( "[BonminOpt::get_variables_types] m != getConstraintCount()" );

    // NOTE: this is hard coded for now, non-lin cases not tested
    for ( size_t i = 0; i != _delegate.getConstraintCount(); ++i )
        //const_types[ i ] = Ipopt::TNLP::LINEAR;
        const_types[ i ] = Ipopt::TNLP::NON_LINEAR;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_constraints_linearity()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_nlp_info( Ipopt::Index                & n
                                   , Ipopt::Index                & m
                                   , Ipopt::Index                & nnz_jac_g
                                   , Ipopt::Index                & nnz_h_lag
                                   , Ipopt::TNLP::IndexStyleEnum & index_style )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ", m = " << m
                  << ")" << std::endl; fflush( stdout );
    }

    n           = _delegate.getVarCount();            // number of variable
    m           = _delegate.getConstraintCount();     // number of constraints
    nnz_jac_g   = _delegate.getJacobian( _ones ).nonZeros(); // number of non zeroes in Jacobian
    VectorX ones( m, 1 ); ones.setConstant( _Scalar(1.) );
    nnz_h_lag   = _delegate.getHessian( 1., ones ).nonZeros();  // number of non zeroes in Hessian of Lagrangean
    index_style = Ipopt::TNLP::C_STYLE;               // zero-indexed

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_nlp_info()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_bounds_info( Ipopt::Index    n
                                      , Ipopt::Number * x_l
                                      , Ipopt::Number * x_u
                                      , Ipopt::Index    m
                                      , Ipopt::Number * g_l
                                      , Ipopt::Number * g_u )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
              << ", m = " << m
              << ")" << std::endl; fflush( stdout );
    }

    if ( n != static_cast<Ipopt::Index>(_delegate.getVarCount()) )
        throw new BonminOptException( "[BonminOpt::get_bounds_info] n != getVarCount()" );
    if ( m != static_cast<Ipopt::Index>(_delegate.getConstraintCount()) )
        throw new BonminOptException( "[BonminOpt::get_bounds_info] m != getConstraintCount()" );

    for ( size_t j = 0; j != _delegate.getVarCount(); ++j )
    {
        x_l[ j ] = _delegate.getVarLowerBound( j );
        x_u[ j ] = _delegate.getVarUpperBound( j );
    } //...for vars

    for ( size_t i = 0; i != _delegate.getConstraintCount(); ++i )
    {
        g_l[ i ] = _delegate.getConstraintLowerBound( i );
        g_u[ i ] = _delegate.getConstraintUpperBound( i );
    } //...for constraints

//    x_l[0] = 0;
//    x_u[0] = 1;

//    x_l[1] = 0;
//    x_u[1] = 1;

//    x_l[2] = 0;
//    x_u[2] = 1;

//    x_l[3] = 0;
//    x_u[3] = 1;

//    g_l[0] = 1.;
//    g_u[0] = DBL_MAX;

//    g_l[1] = 1.;
//    g_u[1] = DBL_MAX;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl; fflush( stdout );
    }

    return true;
} //...BonminOpt::get_bounds_info()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::get_starting_point( Ipopt::Index    n
                                         , bool            init_x
                                         , Ipopt::Number * x
                                         , bool            init_z
                                         , Ipopt::Number * z_L
                                         , Ipopt::Number * z_U
                                         , Ipopt::Index    m
                                         , bool            init_lambda
                                         , Ipopt::Number * lambda )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ", m = " << m
                  << ")" << std::endl; fflush( stdout );
    }

    if ( n != static_cast<Ipopt::Index>(_delegate.getVarCount()) )
        throw new BonminOptException( "[BonminOpt::get_starting_point] n != getVarCount()" );
    if ( m != static_cast<Ipopt::Index>(_delegate.getConstraintCount()) )
        throw new BonminOptException( "[BonminOpt::get_starting_point] m != getConstraintCount()" );

    if ( !init_x )
        throw new BonminOptException( "[BonminOpt::get_starting_point] !init_x..." );
    //assert( init_x );

    if ( init_lambda )
        throw new BonminOptException( "[BonminOpt::get_starting_point] init_lambda..." );
    if ( init_z )
        throw new BonminOptException( "[BonminOpt::get_starting_point] init_z..." );
    //assert( !init_lambda );

    // random for now...
    if ( _delegate.isUseStartingPoint() )
    {
        if ( _delegate.isDebug() )
        {
            std::cout << "using starting point " << _delegate.getStartingPoint().transpose() << std::endl;
        }

        for ( size_t j = 0; j != _delegate.getVarCount(); ++j )
        {
            x[j] = _delegate.getStartingPoint()(j);
        }
    }
    else
    {
        std::random_device  rd;
        std::mt19937        gen(rd());
        std::uniform_int_distribution<>         *int_distribution  = NULL;
        std::uniform_real_distribution<_Scalar> *real_distribution = NULL;
        for ( size_t j = 0; j != _delegate.getVarCount(); ++j )
        {
            if (    (_delegate.getVarType(j) == ParentType::VAR_TYPE::INTEGER)
                 || (_delegate.getVarType(j) == ParentType::VAR_TYPE::BINARY ) )
            {
                if ( !int_distribution )
                    int_distribution = new std::uniform_int_distribution<>( _delegate.getVarLowerBound(j), _delegate.getVarUpperBound(j) );
                x[ j ] = (*int_distribution)( gen );
            }
            else
            {
                if ( !real_distribution )
                {
                    if ( _delegate.isDebug() )
                    {
                        std::cout << "[" << __func__ << "]: " << "creating uniform_real_distributionn with bounds " << _delegate.getVarLowerBound(j) << "..." << _delegate.getVarUpperBound(j) << std::endl;
                    }
                    real_distribution = new std::uniform_real_distribution<_Scalar>( _delegate.getVarLowerBound(j), _delegate.getVarUpperBound(j) );
                }

                x[ j ] = (*real_distribution)( gen );

                if ( std::isinf(x[j]) )
                    std::cerr << "[" << __func__ << "]: " << "warning, returning inf as starting point x0[" << j << "]..." << std::endl;

                if ( _delegate.isDebug() )
                {
                    std::cout << "[" << __func__ << "]: " << "rolled random: " << x[j] << std::endl;
                }
            }
        } //...for variables

        //    x[0] = 1;
        //    x[1] = 1;
        //    x[2] = 1;
        //    x[3] = 1;

        // cleanup
        if ( int_distribution ) { delete int_distribution; int_distribution = NULL; }
        if ( real_distribution ) { delete real_distribution; real_distribution = NULL; }
    }

    if ( _delegate.isDebug() )
    {
        std::cout<<"x:";for(Ipopt::Index vi=0;vi!=static_cast<Ipopt::Index>(n);++vi)std::cout<<x[vi]<<" ";std::cout << "\n";
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
        fflush( stdout );
    }

    return true;
} //...BonminOpt::get_starting_point()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_f( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value )
{
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ")" << std::endl;
        fflush( stdout );
    }

    if ( n != static_cast<Ipopt::Index>(_delegate.getVarCount()) )
        throw new BonminOptException( "[BonminOpt::eval_f] n != getVarCount()" );


    MatrixMapT x_eig( x, _delegate.getVarCount() );
    obj_value = (x_eig.transpose() * _delegate.getCachedQo() * x_eig + x_eig.transpose() * _delegate.getCachedqo() ).coeff( 0 );

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "obj_value: " << obj_value << " from x " << x_eig.transpose() << std::endl;
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
        fflush( stdout );
    }

    return true;
} //...BonminOpt::eval_f()

#define TIC auto start = std::chrono::system_clock::now();
#define RETIC start = std::chrono::system_clock::now();
#define TOC(title,it) { std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start; \
                     std::cout << title << ": " << elapsed_seconds.count()/it << " s" << std::endl; }

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_grad_f( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
#if QCQP_DEBUG
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ")" << std::endl;
        fflush( stdout );
    }
#endif

#if 0
    //if ( new_x || !_delegate.getGradF().rows() )
    {
        if ( n != _delegate.getVarCount() )
            throw new BonminOptException( "[BonminOpt::eval_grad_f] n != getVarCount()" );

        // first derivatives of linear objective function are the coefficients themselves
        //Eigen::Matrix<_Scalar,Eigen::Dynamic,1> _grad_f( _delegate.getCachedqo() );
        _delegate.getGradF() = _delegate.getCachedqo();

        //SparseMatrix const& Qo          ( _delegate.getCachedQo() ); // cached quadratic matrix
        MatrixMapT          x_eig       ( x, n );      // input x

        _delegate.getGradF() += (_delegate.getCachedQo() * x_eig).eval() + (_delegate.getCachedQo().transpose() * x_eig).eval();
    }

    // copy to output
    MatrixNonConstMapT( grad_f, n ) = _delegate.getGradF();
#else
    MatrixNonConstMapT grad_f_eig( grad_f, n );
    MatrixMapT         x_eig     ( x, n );      // input x

    grad_f_eig = _delegate.getCachedqo() +
                 (_delegate.getCachedQo()             * x_eig).eval() +
                 (_delegate.getCachedQo().transpose() * x_eig).eval();
#endif


#if QCQP_DEBUG
  if ( _delegate.isDebug() )
  {
      std::cout << "[" << __func__ << "]: " << "grad_f: ";
      for(size_t vi=0;vi!=n;++vi)std::cout<<grad_f[vi]<<" ";std::cout << "\n"; std::cout<<"x:";for(size_t vi=0;vi!=n;++vi)std::cout<<x[vi]<<" ";std::cout << "\n";
      std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
      fflush( stdout );
  }
#endif

  return true;
} //...BonminOpt::eval_grad_f()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_g( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g )
{

#   if QCQP_DEBUG
//    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ", m = " << m
                  << ")" << std::endl;
        fflush( stdout );
    }
#   endif

    if ( n != static_cast<Ipopt::Index>(_delegate.getVarCount()) )
        throw new BonminOptException( "[BonminOpt::eval_g] n != getVarCount()" );
    if ( m != static_cast<Ipopt::Index>(_delegate.getConstraintCount()) )
        throw new BonminOptException( "[BonminOpt::eval_g] m != getConstraintCount()" );

    MatrixMapT                  x_eig( x, _delegate.getVarCount() );
    Eigen::Matrix<_Scalar,-1,1> c    ( _delegate.getCachedA() * x_eig );

    for ( size_t j = 0; j != _delegate.getConstraintCount(); ++j )
    {
        if ( (j < _delegate.getCachedQSize()) && (_delegate.getCachedQ(j).nonZeros()) )
            c(j) += (x_eig.transpose() * _delegate.getCachedQ(j) * x_eig).coeff(0);

        g[j] = c(j);
    }

#   if QCQP_DEBUG
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
        fflush( stdout );
    }
#   endif

    return true;
} //...BonminOpt::eval_g()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_jac_g( Ipopt::Index         n
                              , Ipopt::Number const* x
                              , bool          new_x
                              , Ipopt::Index         m
                              , Ipopt::Index         nnz_jac
                              , Ipopt::Index       * iRow
                              , Ipopt::Index       * jCol
                              , Ipopt::Number      * values )
{
    bool ret_val = false;

#if QCQP_DEBUG
    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
                  << ", m = " << m
                  << ")" << std::endl;
        fflush( stdout );
    }
#endif

    if ( n != static_cast<Ipopt::Index>(_delegate.getVarCount()) )
    {
        std::cerr  << "[" << __func__ << "]: " << "n " << n << " != " << _delegate.getVarCount() << " getVarCount()\n";
        throw new BonminOptException( "[BonminOpt::eval_jac_g] n != getVarCount()" );
    }

    MatrixMapT   x_eig( x, _delegate.getVarCount() );
    SparseMatrix jacobian = x ? _delegate.getJacobian( x_eig )
                              : _delegate.getJacobian( _ones );

    if ( nnz_jac != jacobian.nonZeros() )
        throw new BonminOptException( "[BonminOpt::eval_jac_g] nnz_jac != _jacobian.nonZeros()" );

    int entry_id = 0;
    if ( values == NULL )
    {
        for ( int row = 0; row != jacobian.outerSize(); ++row )
        {
            for ( typename Eigen::SparseMatrix<_Scalar,Eigen::RowMajor>::InnerIterator it(jacobian,row);
                  it; ++it, ++entry_id )
            {
                iRow[ entry_id ] = it.row();
                jCol[ entry_id ] = it.col();
            } // ...for col
        } // ...for row

        ret_val = true;
    } // ... if values == NULL
    else
    {
        for ( int row = 0; row != jacobian.outerSize(); ++row )
        {
            for ( typename Eigen::SparseMatrix<_Scalar,Eigen::RowMajor>::InnerIterator it(jacobian,row);
                  it; ++it, ++entry_id )
            {
                values[ entry_id ] = it.value();
            } // ...for col
        } // ...for row

        ret_val = true;
    } // ... else values != NULL

#if QCQP_DEBUG
    if ( _delegate.isDebug() )
    {
        //std::cout << "[" << __func__ << "]: " << "returned " << jacobian;
        //if ( x ) std::cout << "\nfrom " << MatrixMapT(x, _delegate.getVarCount()).transpose();
        //std::cout << std::endl;

        std::cout << "[" << __func__ << "]: " << "finish" << std::endl;
        fflush( stdout );
    }
#endif

    return ret_val;
} // ... BonminOpt::eval_jac_g()

template <typename _Scalar> bool
BonminTMINLP<_Scalar>::eval_h( Ipopt::Index          n
                             , Ipopt::Number const * x
                             , bool           new_x
                             , Ipopt::Number         obj_factor
                             , Ipopt::Index          m
                             , Ipopt::Number const * lambda
                             , bool           new_lambda
                             , Ipopt::Index          nele_hess
                             , Ipopt::Index        * iRow
                             , Ipopt::Index        * jCol
                             , Ipopt::Number       * values )
{
    bool ret_val = false;

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "call(n = " << n
              << ", m = " << m
              << ")" << std::endl;
        fflush( stdout );
    }


    if ( n != static_cast<Ipopt::Index>(_delegate.getVarCount()) )
        throw new BonminOptException( "[BonminOpt::eval_h] n != getVarCount()" );
    if ( m != static_cast<Ipopt::Index>(_delegate.getConstraintCount()) )
        throw new BonminOptException( "[BonminOpt::eval_h] m != getConstraintCount()" );


    SparseMatrix hessian = lambda ? _delegate.getHessian( obj_factor, MatrixMapT(lambda, m) )
                                  : _delegate.getHessian( obj_factor, VectorX::Ones(m,1) );

    if ( nele_hess != hessian.nonZeros() )
        throw new BonminOptException( "[BonminOpt::eval_h] nele_hess != _delegate.getHessian().nonZeros()" );

    int entry_id = 0;
    if ( values == NULL )
    {
        for ( int row = 0; row != hessian.outerSize(); ++row )
        {
            for ( typename Eigen::SparseMatrix<_Scalar,Eigen::RowMajor>::InnerIterator it(hessian,row);
                  it; ++it, ++entry_id )
            {
                iRow[ entry_id ] = it.row();
                jCol[ entry_id ] = it.col();
            } // ...for col
        } // ...for row

        ret_val = true;
    }
    else {
        // NOTE: lower triangular only please!
        for ( int row = 0; row != hessian.outerSize(); ++row )
        {
            for ( typename Eigen::SparseMatrix<_Scalar,Eigen::RowMajor>::InnerIterator it(hessian,row);
                  it; ++it, ++entry_id )
            {
                values[ entry_id ] = it.value() * obj_factor;
            } // ...for col
        } // ...for row
        ret_val = true;
    }

    if ( _delegate.isDebug() )
    {
        std::cout << "[" << __func__ << "]: " << "finished" << std::endl;
        fflush( stdout );
    }

    return ret_val;
} // ...BonminOpt::eval_h()

template <typename _Scalar> void
BonminTMINLP<_Scalar>::finalize_solution( TMINLP::SolverReturn   status
                                        , Ipopt::Index                  n
                                        , Ipopt::Number const         * x
                                        , Ipopt::Number                 obj_value )
{
    std::cout << "Problem status: "  << status;
    switch (status)
    {
        case SUCCESS: std::cout << " SUCCESS"; break;
        case INFEASIBLE: std::cout << " INFEASIBLE"; break;
        case CONTINUOUS_UNBOUNDED: std::cout << " CONTINUOUS_UNBOUNDED"; break;
        case LIMIT_EXCEEDED: std::cout << " LIMIT_EXCEEDED"; break;
        case MINLP_ERROR: std::cout << " MINLP_ERROR"; break;
        default: std::cout << " UNKNOWN"; break;
    }
    std::cout << std::endl;

    std::cout << "Objective value: " << obj_value << std::endl;
    if ( _delegate.isPrintSol() && (x != NULL) )
    {
        std::cout << "Solution:" << std::endl;
        for ( int i = 0 ; i < n; ++i )
        {
            std::cout << "x[" << i << "] = " << x[i];
            if ( i < n-1 ) std::cout << ", ";
        }
        std::cout << std::endl;
    } // ... printSol_

    // save solution
    VectorX sol( n );
    sol.setZero();
    if ( x )
    {
        std::cout << "saving " << n << " variables" << std::endl;
        _Scalar sum = 0;
        for ( int i = 0 ; i < n; ++i )
        {
            sol(i) = x[i];
            sum += x[i];
        }
        std::cout << "sum: " << sum << std::endl;
    }
    else
    {
        std::cerr << "[" << __func__ << "]: " << "no solution returned, status is ";
        switch ( status )
        {
            case    SUCCESS: std::cerr << "SUCCESS"; break;
            case INFEASIBLE: std::cerr << "INFEASIBLE"; break;
            case CONTINUOUS_UNBOUNDED: std::cerr << "CONTINUOUS_UNBOUNDED"; break;
            case LIMIT_EXCEEDED: std::cerr << "LIMIT_EXCEEDED"; break;
            case MINLP_ERROR: std::cerr << "MINLP_ERROR"; break;
            default: std::cerr << "unrecognized code"; break;
        }
        std::cerr << std::endl;
        fflush( stderr );

    } //...if x

    _delegate.setSolution( sol );
} // ... BonminOpt::finalize_solution()

} //...ns GF2

#endif // BONMINOPT_H
