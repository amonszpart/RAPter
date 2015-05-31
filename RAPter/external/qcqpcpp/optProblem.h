#ifndef QCQPCPP_INC_SGOPTPROBLEM_H
#define QCQPCPP_INC_SGOPTPROBLEM_H

#include <iostream> // cout, cerr, endl
#include "Eigen/Sparse"

namespace qcqpcpp
{

class OptProblemException : public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

template <typename _Scalar>
class OptProblem
{
    public:
        typedef _Scalar                                     Scalar;         //!< \brief Some implementations require this to be double (i.e. Mosek does).
        typedef int                                         ReturnType;     //!< \brief Error code storage type.
        typedef Eigen::Matrix<_Scalar,-1,1>                 VectorX;        //!< \brief ColVector of unkown dimensions.
        typedef Eigen::Triplet<Scalar>                      SparseEntry;    //!< \brief This type is an element of a list of entries before construction of a SparseMatrix.
        typedef std::vector<SparseEntry>                    SparseEntries;  //!< \brief A list of entries before construction. SparseMatrices are created in a lazy fashion.
        typedef Eigen::SparseVector<Scalar,Eigen::RowMajor> SparseVector;   //!< \brief Unused.
        typedef Eigen::SparseMatrix<Scalar,Eigen::RowMajor> SparseMatrix;   //!< \brief To store quadratic objective (Q_o) and quadratic constraints (Qi). SparseMatrices are created in a lazy fashion.

        enum OBJ_SENSE   { MINIMIZE       //!< \brief Minimize objective function in #optimize()
                         , MAXIMIZE       //!< \brief Maximize objective function in #optimize()
                         };
        enum VAR_TYPE    { CONTINUOUS = 0 //!< \brief Default
                         , INTEGER    = 1 //!< \brief Mixed Integer Programming
                         , BINARY     = 2 //!< \brief Binary integer programming. Some solvers allow this option to be set.
                         };
        enum LINEARITY   { LINEAR     = 0 //!< \brief Default value.
                         , NON_LINEAR = 1 //!< \brief NOT USED, here for completeness sake
                         };
        enum BOUND       { GREATER_EQ = 0 //!< \brief Means: lower_bound < ... < +INF
                         , LESS_EQ    = 1 //!< \brief Means: -INF        < ... < upper_bound
                         , EQUAL      = 2 //!< \brief Means: Const       = ... = Const
                         , FREE       = 3 //!< \brief Means: -INF        = ... = +INF
                         , RANGE      = 4 //!< \brief Means: lower_bound = ... = upper_bound
                         };

        //! \brief Internal structure to store a variable
        struct Variable
        {
                BOUND          _bkx;                 //!< \brief Variable bound key
                Scalar         _blx;                 //!< \brief Variable numerical value of lower bound
                Scalar         _bux;                 //!< \brief Variable numerical value of upper bound
                VAR_TYPE       _type_x;              //!< \brief Variable type
                LINEARITY      _lin_x;               //!< \brief Vairable linearity.
                std::string    _name;                //!< \brief Name of variable for human reasons.
        };

        //! \brief Virtual getter for the implementation specific infinity value. Used, when a bound direction is "unbounded"
        //! \return The value infinity, i.e. DBL_MAX or MSK_INFINITY.
        virtual Scalar                           getINF                 ()                const { return std::numeric_limits<Scalar>::max(); }
        //!< \brief Virtual getter for implementation specific error code for "no error".
        //! \return The int conversion of the code (i.e. 0).
        virtual int                              getOkCode              ()                const { return 0; }

        //! \brief Make sure, to call before optimize is called, but AFTER the problem details are added using the setters below.
        virtual ReturnType                       update                 ( bool verbose = false )                                          { std::cerr << "[" << __func__ << "]: " << "Please override with library specific set-up logic. See MosekOpt.h for concept." << std::endl; return ReturnType(EXIT_FAILURE); }
        //! \brief Call to run optimization, AFTER the problem details were added, and <i>update()</i> was called.
        virtual ReturnType                       optimize               ( std::vector<Scalar> *x_out           = NULL,
                                                                          OBJ_SENSE            objecitve_sense = OBJ_SENSE::MINIMIZE ) { std::cerr << "[" << __func__ << "]: " << "Please override with library specific set-up logic. See MosekOpt.h for concept." << std::endl; return ReturnType(EXIT_FAILURE); }

        //! \brief Constructor unused for the moment. Start with <i>addVariable()</i>.
                                                  OptProblem             () : _updated(false), _useStartingPoint(false), _time_limit(Scalar(-1)), _tol_rel_gap(Scalar(-1)) {}
        //! \brief Destructor unused for the moment. Declared virtual for inheritence.
        virtual                                  ~OptProblem             () { /*std::cout << "[" << __func__ << "][INFO]: " << "Empty destructor" << std::endl;*/ }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Variables ////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //! \brief addVariable  Start setting up the problem by defining the X vector.
        //! \param bound_type   A variable can be FREE (unbounded), EQUAL (fixed), lower bounded (GREATER_EQ) or upper bounded (LESS_EQ), or bounded (RANGE).
        //! \param lower_bound  Specify lower bound (set to -INF), if not lower bounded (FREE, LESS_EQ).
        //! \param upper_bound  Specify upper bound (set to +INF), if not upper bounded (FREE, GREATER_EQ).
        //! \param var_type     Continuous (default) or integer (for mixed integer programming).
        //! \return             0-indexed id of the just added variable.
        inline int                               addVariable            ( BOUND       bound_type     = BOUND::FREE
                                                                        , Scalar      lower_bound    = -getINF()
                                                                        , Scalar      upper_bound    = +getINF()
                                                                        , VAR_TYPE    var_type       = VAR_TYPE::CONTINUOUS
                                                                        , LINEARITY   var_lin        = LINEARITY::LINEAR
                                                                        , std::string name           = "" );
        inline int                               addVariable            ( Variable const& var ) { return this->addVariable( var._bkx, var._blx, var._bux, var._type_x, var._lin_x, var._name ); }

        inline VAR_TYPE                   const  getVarType             ( int j )         const { return _type_x[j]; }
        inline BOUND                      const  getVarBoundType        ( int j )         const { return _bkx[j]; }
        inline Scalar                     const  getVarLowerBound       ( int j )         const { return _blx[j]; }
        inline Scalar                     const  getVarUpperBound       ( int j )         const { return _bux[j]; }
        inline LINEARITY                  const  getVarLinearity        ( int j )         const { return _lin_x[j]; }
        inline size_t                            getVarCount            ()                const { return _bkx.size(); }     //!< \brief Returns the number of variables currently in the system
        inline std::string                       getVarName             ( int var_id )    const { return _names[var_id]; }  //!< \brief Returns variable name.

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Objective function ///////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /// Formulation: minimize/maximize X' * 1/2 * Q_o * X + q_o * X + c

        // Const objective function (c)
        inline int                               setObjectiveBias               ( Scalar cfix ) { _cfix = cfix; return EXIT_SUCCESS; } //!< \brief Is always added to the objective function as bias.
        inline Scalar                            getObjectiveBias               ()        const { return _cfix; }                      //!< \brief Read-only getter for fixed objective function bias.
        // Linear objective function (q_o)
        inline int                               setLinObjective                ( int j, Scalar coeff );                 //!< \brief Sets the j-th linear coefficient to coeff, regardless previous value.
        inline int                               setLinObjectives               ( Eigen::Matrix<_Scalar,-1,1> const& obj ) { for ( int i = 0; i != obj.rows(); ++i ) setLinObjective(i, obj(i) ); return EXIT_SUCCESS; }
        inline int                               addLinObjective                ( int j, Scalar coeff );                 //!< \brief Adds coeff to the j-th linear coefficient.
        inline int                               addLinObjectives               ( SparseMatrix const& mx );              //!< \brief Assumes vector input, so ((mx.cols == 1) || (mx.rows == 1)).
        inline std::vector<Scalar>        const& getLinObjectives               ()        const { return _linObjs; }     //!< \brief Read-only getter for linear objective vector.
        inline SparseMatrix                      getLinObjectivesMatrix         ()        const;                         //!< \brief Returns a 1 column SparseMatrix temporary.
        // Quadratic Objecitve function (Q_o)
        inline int                               addQObjective                  ( size_t i, size_t j, Scalar coeff );          //!< \brief Adds input to quadratic objective matrix, duplicate entries get summed in a lazy fashion (when update() is called).
        inline int                               addQObjectives                 ( SparseMatrix const& mx );              //!< \brief Adds input matrix to quadratic objective matrix, duplicate entries get summed in a lazy fashion (when update() is called).
        inline int                               setQObjectives                 ( SparseMatrix const& mx );              //!< \brief Clears current objectives, and sets input as quadratic objective matrix.
        inline SparseEntries              const& getQuadraticObjectives         ()        const { return _quadObjList; } //!< \brief Read-only getter for quadratic objective matrix's unordered entry set.
        inline SparseMatrix                      getQuadraticObjectivesMatrix   ()        const;                         //!< \brief Returns all Qo entries assembled to a Sparsematrix, temporarily.
        // Hessian of Lagrangian
        //! \brief  estimateHessianOfLagrangian  Estimates Hessian of Lagrangian from Quadratic objective matrix.
        //! \return SparseMatrix                 Lower triangle only, assumed to be symmetric. (TODO: add quadratic constraints).
        //inline SparseMatrix                      estimateHessianOfLagrangian    ()        const;
        inline int                               precalcHessianCoeffs           ();
        inline std::vector<SparseMatrix>  const& getHessians                    ()        const { return _hessians; }
        inline SparseMatrix               const  getHessian                     ( Scalar const obj_factor, VectorX const& lambdas ) const;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Constraints //////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /// Formulation: lower_bound_i <= 1/2 * X' * Qi * X + A(i,:) * X <= upper_bound_i
        /// Note: There can only be one objective function, but there can be m constraints, indexed by i.
        /// TODO: assumed, not tested yet...
        /// Tested:   Q0     AND A with one row (A vector)
        /// Tested:              A matrix alone
        /// UNtested: Q0, Q1 AND A matrix

        //! \brief                      Append a constraint to the system, and optionally a linear line to the constraint matrix A. The line's row_id is the same, as the constr_id of the quadratic constraints.
        //! \param bound_type           The bound type of the i-th constraint (i-th constraint is assumed to mean: Qi and A(i,:) )
        //! \param lower_bound          The lower bound of the i-th constraint (i-th constraint is assumed to mean: Qi and A(i,:) )
        //! \param upper_bound          The upper bound of the i-th constraint (i-th constraint is assumed to mean: Qi and A(i,:) )
        //! \param coeffs               Pointer to a dense vector to be appended to linear constraint matrix A.
        //! \return                     EXIT_SUCCESS/EXIT_FAILURE, if bound!=SG_INFINITY was set to an unbounded type (i.e. SG_BOUND::FREE, lower bound of SG_BOUND::LESS_EQ, or upper bound of SG_BOUND::GREATER_EQ)
        inline int                               addConstraint                  ( BOUND                bound_type
                                                                                , Scalar               lower_bound // = -getINF()
                                                                                , Scalar               upper_bound // = +getINF()
                                                                                , SparseMatrix       * row_vector     = NULL
                                                                                , LINEARITY            c_lin = LINEARITY::LINEAR
                                                                                );
        inline int                               addConstraint                  ( Variable const& constr );

        //! \brief                      Add linear constraint using an std::vector
        inline int                               addLinConstraint               ( BOUND                bound_type
                                                                                , Scalar               lower_bound // = -getINF()
                                                                                , Scalar               upper_bound //  = +getINF()
                                                                                , std::vector<Scalar>  coeffs
                                                                                , LINEARITY            c_lin = LINEARITY::LINEAR
                                                                                );

        inline int                               addLinConstraints              ( SparseMatrix const& mx );                     //!< \brief Adds full linear constraint matrix. Please call addLinConstraint before.
        inline int                               setLinConstraints              ( SparseMatrix const& mx );                     //!< \brief Sets full linear constraint matrix, clears previous lin constraints. Please call addLinConstraint before.

        inline BOUND                      const  getConstraintBoundType         ( int i ) const { return _bkc[i]; }             //!< \brief Returns bound type of i-th constraint (free, lower, upper, range, fixed).
        inline Scalar                     const  getConstraintLowerBound        ( int i ) const { return _blc[i]; }             //!< \brief Returns lower bound of i-th constraint.
        inline Scalar                     const  getConstraintUpperBound        ( int i ) const { return _buc[i]; }             //!< \brief Returns upper bound of i-th constraint.
        inline LINEARITY                  const  getConstraintLinearity         ( int j ) const { return _lin_c[j]; }
        inline size_t                            getConstraintCount             ()        const { return _bkc.size(); }         //!< \brief Returns the number of linear constraint lines currently in the system.
        //inline SparseMatrix                      estimateJacobianOfConstraints  ()        const;
        inline SparseMatrix               const  getJacobian                    ( VectorX const& x );                           //!< \brief Calculated on the fly from input vector x and precalculated _jacobian.

        //! \brief addQConstraint       Append an entry to a quadratic constraint matrix. Duplicate entries are summed in a lazy fashion (when <i>update()</i>) is called.
        //! \brief                      The corresponding linear constraints and bounds are stored in A.row( constr_id ), and _b[k|l|u]c[constr_id].
        /// \param constr_id            Append an entry to Qi, where i == constr_id. Create entry by <i>addLinConstraint</i> before a quadratic entry can be added.
        /// \param i                    Row index of entry in Qi.
        /// \param j                    Column index of entry in Qi.
        /// \param coeff                Entry value to insert into Qi.
        /// \return                     EXIT_SUCCESS/EXIT_FAILURE, if constr_id >= getConstraintCount().
        inline int                               addQConstraint                 ( size_t constr_id, size_t i, size_t j, Scalar coeff );
        inline int                               addQConstraints                ( SparseMatrix const& mx );

        inline SparseEntries              const& getLinConstraints              ()        const { return _linConstrList; }      //!< \brief Returns all A entries, unordered.
        inline SparseMatrix                      getLinConstraintsMatrix        ()        const;                                //!< \brief Returns all A entries assembled to a Sparsematrix, temporarily.
        inline SparseEntries              const& getQuadraticConstraints        ( int i ) const { return _quadConstrList[i]; }  //!< \brief Returns all Qi entries, unordered.
        inline std::vector<SparseEntries> const& getQuadraticConstraints        ()        const { return _quadConstrList; }     //!< \brief Returns vector of Qi entries, vector is ordered by i, but the entries are unordered.
        inline SparseMatrix                      getQuadraticConstraintsMatrix  ( size_t i ) const;                             //!< \brief Returns all Qi entries assembled to a Sparsematrix, temporarily.

        //! \brief              Prints inner state before optimization.
        //! \param entry_limit  How many entries to print of each matrix.
        inline int                               printProblem                   ( size_t entry_limit = 5 ) const;

        //! \brief Setter to store latest solution, to be called from #optimize().
        inline void                              setSolution                    ( VectorX const& sol ) { _x = sol; }

        //! \brief Setter to set initial solution.
        inline int                               setStartingPointDense          ( VectorX      const& x0 ) { _x0 = x0; _useStartingPoint = true; return EXIT_SUCCESS; }
        inline int                               setStartingPoint               ( SparseMatrix const& x0 );
        //! \brief Getter to check, if _x0 should be used as starting point for the optimization.
        inline bool                              isUseStartingPoint             ()        const        { return _useStartingPoint; }
        inline VectorX const&                    getStartingPoint               ()        const        { return _x0; }
        SparseMatrix                             getStartingPointMatrix         () const;                                       //!< \brief Returns X0, calculated on the fly.

        // IO
        inline int                               write                          ( std::string const& path ) const;              //!< \brief Serialize to disk.
        inline int                               read                           ( std::string        proj_path );               //!< \brief Read from disk.

        // PARAMS
        //! \brief            Set time limit to allow the implementation to run for.
        //! \param time_limit The time limit to be set. Has to be > 0.
        inline int                               setTimeLimit                   ( Scalar time_limit ) { _time_limit = time_limit; return EXIT_SUCCESS; }
        inline Scalar                            getTimeLimit                   () const              { return _time_limit; }   //!< \brief Returns solver time limit parameter.
        //! \brief             Set relative gap tolerance for the solver.
        //! \param tol_rel_gap The gap tolerance to set.
        inline int                               setTolRelGap                   ( Scalar tol_rel_gap) { _tol_rel_gap = tol_rel_gap; return EXIT_SUCCESS; }
        inline Scalar                            getTolRelGap                   () const              { return _tol_rel_gap; }   //!< \brief Returns solver relative gap tolerance.

        inline std::string                      getAuxName()                              const { return "aux.csv"; } //!< \brief File name constexpr to store variables and constraints. Used by #write().
        inline std::string                      getqoName()                               const { return "qo.csv";  } //!< \brief File name constexpr to store linear objective vector. Used by #write().
        inline std::string                      getQoName()                               const { return "Qo.csv";  } //!< \brief File name constexpr to store quadratic objective matrix. Used by #write().
        inline std::string                      getAName()                                const { return "A.csv";   } //!< \brief File name constexpr to store linear constraints matrix. Used by #write().
        inline std::string                      getX0Name()                               const { return "X0.csv";  } //!< \brief File name constexpr to store starting point (initial solution) vector. Used by #write().
        //! \brief File name constexpr to store qaudratic constraints matrices. Used by #write().
        inline std::string                      getQiName( int i )                        const { char id[255]; sprintf( id, "Q%d.csv", i ); return std::string( id ); }

    protected:
        // Variables: _blx <= X <= _buc
        std::vector<BOUND>          _bkx;                 //!< \brief Variable bound keys
        std::vector<Scalar>         _blx;                 //!< \brief Variable numerical value of lower bounds
        std::vector<Scalar>         _bux;                 //!< \brief Variable numerical value of upper bounds
        std::vector<VAR_TYPE>       _type_x;              //!< \brief Variable type
        std::vector<LINEARITY>      _lin_x;               //!< \brief Variable linearity
        std::vector<std::string>    _names;               //!< \brief Variable names.

        // Objective: 1/2 * X' * Q_o * X + q_o * X + c
        Scalar                      _cfix;                //!< \brief Constant bias to objective function                                   (c)
        std::vector<Scalar>         _linObjs;             //!< \brief Linear variable coeffs (objective function)                           (q_o)
        SparseEntries               _quadObjList;         //!< \brief Accumulator for QObjective values                                     (Q_o)

        // Constraints of the form: _blc <= 1/2 * X' * Qi * X + A(i,:) * X <= _buc //TODO: is it really A(i,:)?
        std::vector<BOUND>          _bkc;                 //!< \brief Constraints bound keys
        std::vector<Scalar>         _blc;                 //!< \brief Constraints numerical value of lower bounds
        std::vector<Scalar>         _buc;                 //!< \brief Constraints numerical value of upper bounds
        std::vector<LINEARITY>      _lin_c;               //!< \brief Constraints linearity
        SparseEntries               _linConstrList;       //!< \brief Stores values to construct the A linear constraint matrix, unordered. (A, A.cols = n = # of vars, A.rows = m = # of constraints)
        std::vector<SparseEntries>  _quadConstrList;      //!< \brief Accumulator for Quadratic constraint values                           [Q0, Q1, ... Qi ... Qm ]
        std::vector<SparseMatrix>   _hessians;            //!< \brief [0] = Hessian of ObjectiveFunction, [1..constrCount] = Hessian of Constraints
        SparseMatrix                _jacobian;            //!< \brief Fixed part (linear) of jacobian of constraints.

        bool                        _updated;             //!< \brief This has to be flipped to true by calling update() for optimize() to run.
        VectorX                     _x;                   //!< \brief Optimize() stores the solution here when finishing.
        VectorX                     _x0;                  //!< \brief Starting point for optimization. TODO: store sparse instead...
        bool                        _useStartingPoint;    //!< \brief Decides, if the starting point should be used (if false, a random start is used.) \warning Unused in Mosek, Gurobi

        // Parameters
        Scalar                      _time_limit;          //!< \brief Time limit to allow implementation to run for.
        Scalar                      _tol_rel_gap;         //!< \brief Relative Gap tolerance.

        inline int                              _parseAuxFile( std::string const& aux_file_path );                    //!< \brief Method to read variables and constraints from file. Used by #read().
}; // ...class SGOpt

} // ... namespace qcqpcpp

#ifndef QCQPCPP_INC_SGOPTPROBLEM_HPP
#   define QCQPCPP_INC_SGOPTPROBLEM_HPP
#   include "qcqpcpp/impl/optProblem.hpp"
#endif // QCQPCPP_INC_SGOPTPROBLEM_HPP

#endif // QCQPCPP_INC_SGOPTPROBLEM_H
