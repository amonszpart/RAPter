#ifndef QCQPCPP_SGOPTPROBLEM_HPP
#define QCQPCPP_SGOPTPROBLEM_HPP

#include <exception>
#include "qcqpcpp/io/io.h" // readSparseMatrix, writeSparseMatrix
#include "sys/stat.h"      // mkdir

namespace qcqpcpp
{

template <class _Derived>
inline void print( std::string const& title, _Derived const& mx )
{
    std::cout << title << std::endl;
    for ( int y = 0; y != mx.rows(); ++y )
    {
        for ( int x = 0; x != mx.cols(); ++x )
        {
            std::cout << mx.coeff(y,x) << " ";
        }
        std::cout << std::endl;
    }
}

template <typename _Scalar> int
OptProblem<_Scalar>::addVariable( BOUND bound_type, Scalar lower_bound, Scalar upper_bound, VAR_TYPE var_type, LINEARITY var_linearity /* = LINEAR */, std::string name /* = "" */ )
{
    // var bound type { MSK_BK_FX = fixed, MSK_BK_FR = free, MSK_BK_LO = blx[j] .. INF, MSK_BK_RA = blx[j] .. bux[j], MSK_BK_UP = -INF .. bux[j] }
    _bkx.push_back( bound_type  );
    // lower bounds
    _blx.push_back( lower_bound );
    // upper bounds
    _bux.push_back( upper_bound );
    // linear objective coeff
    _linObjs  .push_back( Scalar(0)   );
    // variable type
    _type_x.push_back( var_type );
    // variable linearity
    _lin_x.push_back( var_linearity );
    // variable name
    _names.push_back( name );

    //return EXIT_SUCCESS;
    return _bkx.size()-1;
} // ...OptProblem::addVariable

//______________________________________________________________________________

template <typename _Scalar> int
OptProblem<_Scalar>::setLinObjective( int j, Scalar coeff )
{
    if ( static_cast<int>(_linObjs.size()) <= j )
    {
        std::cerr << "[" << __func__ << "]: " << "please add var " << j << " before setting any coeffs" << std::endl;
        return EXIT_FAILURE;
    }

    _linObjs[ j ] = coeff;

    return EXIT_SUCCESS;
} // ...OptProblem::setVarLinCoeff

template <typename _Scalar> int
OptProblem<_Scalar>::addLinObjective( int j, Scalar coeff )
{
    if ( static_cast<int>(_linObjs.size()) <= j )
    {
        std::cerr << "[" << __func__ << "]: " << "please add var " << j << " before adding any coeffs" << std::endl;
        return EXIT_FAILURE;
    }

    _linObjs[ j ] += coeff;

    return EXIT_SUCCESS;
} // ...OptProblem::addVarLinCoeff

template <typename _Scalar> int
OptProblem<_Scalar>::addLinObjectives( SparseMatrix const& mx )
{
    typedef typename SparseMatrix::Index Index;

    int err = EXIT_SUCCESS;

    // mx should be a vector, one dim should be varCount long, the other 1 long.
    if (  ( (mx.outerSize() != static_cast<Index>(getVarCount())) && (mx.innerSize() != static_cast<Index>(getVarCount())) )
       || ( (mx.outerSize() != static_cast<Index>(1            )) && (mx.innerSize() != static_cast<Index>(1            )) )
       )
    {
        std::cerr << "[" << __func__ << "]: " << "getVarCount " << getVarCount() << " != " << std::max(mx.outerSize(),mx.innerSize()) << " variables in mx" << std::endl;
        return EXIT_FAILURE;
    }

    // set quadratic objective
    for ( int k = 0; k < mx.outerSize(); ++k )
        for ( typename OptProblem<_Scalar>::SparseMatrix::InnerIterator it(mx, k); it; ++it )
        {
            //std::cout << "[" << __func__ << "]: " << "adding at " << std::max(it.row(),it.col()) << " = " << it.value() << std::endl;
            err += addLinObjective( std::max(it.row(),it.col()), it.value() );
        }

    return err;
} // ...OptProblem::addVarLinCoeff

template <typename _Scalar> typename OptProblem<_Scalar>::SparseMatrix
OptProblem<_Scalar>::getLinObjectivesMatrix() const
{
    SparseMatrix smx( this->getVarCount(), 1 );
    for ( size_t i = 0; i != this->_linObjs.size(); ++i )
    {
        smx.insert( i, 0 ) = this->_linObjs[i];
    }

    return smx;
} // ...OptProblem::getQuadraticObjectivesMatrix()

template <typename _Scalar> int
OptProblem<_Scalar>::addQObjective( size_t i, size_t j, Scalar coeff )
{
    if ( (i >= this->getVarCount()) || j >= (this->getVarCount()) )
    {
        std::cerr << "[" << __func__ << "]: " << "i " << i << " or j " << j << " > " << this->getVarCount() << ", please call addVariable() first! returning..." << std::endl;
        return EXIT_FAILURE;
    }

    // only lower triangle
    if ( j > i )
        std::swap( i, j );

    _quadObjList.push_back( Eigen::Triplet<Scalar>(i,j,coeff) );

    return EXIT_SUCCESS;
} // ...OptProblem::addQObjectives()

template <typename _Scalar> int
OptProblem<_Scalar>::addQObjectives( SparseMatrix const& mx )
{
    int err = EXIT_SUCCESS;

    typename SparseMatrix::Index varCount = getVarCount();
    if ( (mx.outerSize() != varCount) || (mx.innerSize() != varCount) )
    {
        std::cerr << "[" << __func__ << "]: " << "getVarCount " << getVarCount() << " != " << mx.outerSize() << " || " << mx.innerSize() << " variables in mx" << std::endl;
        return EXIT_FAILURE;
    }

    // set quadratic objective
    _quadObjList.reserve( _quadObjList.size() + mx.nonZeros() );
    for ( int k = 0; k < mx.outerSize(); ++k )
        for ( typename OptProblem<_Scalar>::SparseMatrix::InnerIterator it(mx, k); it; ++it )
            err += addQObjective( it.row(), it.col(), it.value() );

    return err;
} // ...OptProblem::addQObjectives()

template <typename _Scalar> int
OptProblem<_Scalar>::setQObjectives( SparseMatrix const& mx )
{
    _quadObjList.clear();
    return addQObjectives( mx );
} // ...OptProblem::setQObjectives()

template <typename _Scalar> typename OptProblem<_Scalar>::SparseMatrix
OptProblem<_Scalar>::getQuadraticObjectivesMatrix() const
{
    SparseMatrix mx( this->getVarCount(), this->getVarCount() );
    mx.setFromTriplets( this->_quadObjList.begin(), this->_quadObjList.end() );

    return mx;
} // ...OptProblem::getQuadraticObjectivesMatrix()

template <typename _Scalar> typename OptProblem<_Scalar>::SparseMatrix const
OptProblem<_Scalar>::getJacobian( VectorX const& x )
{
    //typedef typename BonminOpt<_Scalar>::SparseEntry SparseEntry;

    // init linear cache part
    if ( !this->_jacobian.size() )
    {
        this->_jacobian = SparseMatrix( getConstraintCount(), getVarCount() );
        this->_jacobian.setFromTriplets( _linConstrList.begin(), _linConstrList.end() );
    }

    // linear part is cached
    SparseMatrix jacobian = this->_jacobian;
    if ( x.rows() != jacobian.cols() )
    {
        std::cerr << "[" << __func__ << "]: "
                  << "input vector x is not varCount sized " << x.rows() << " x " << x.cols() << " .rows != " << jacobian.rows() << " x " << jacobian.cols() << ".cols" << std::endl;
        throw new OptProblemException( "input vector x is not varCount sized" );
    }

    // quadratic part is estimated on the fly
    // plus the coeff * other term from its quadratic matrix
    for ( size_t constr_id = 0; constr_id != this->_quadConstrList.size(); ++constr_id )
        for ( size_t entry_id = 0; entry_id != this->_quadConstrList[constr_id].size(); ++entry_id )
        {
            SparseEntry const& entry = this->_quadConstrList[constr_id][entry_id];

            if ( entry.row() < entry.col() ) // above diagonal
            {
                std::cerr << "[" << __func__ << "]: " << "non-lower-triangular constraint matrix..." << std::endl;
                throw new OptProblemException("non-lower-triangular constraint matrix...");
            }
            else if ( entry.row() == entry.col() ) // on diagonal
            {
                jacobian.coeffRef( /* row == constr_id: */ constr_id
                                 , /* col == var_id:    */ entry.row()                // J(i,j) =
                                 ) += entry.value() * _Scalar(2.) * x( entry.col() ); //        = 2 * c_jj * x_j
            }
            else // below diagonal -> derive to two columns for the two variables
            {
                // d(c_21 * x_2 * x_1) / dx_2 = c_21 * x_1
                jacobian.coeffRef( /* row == constr_id: */ constr_id
                                 , /* col == var_id:    */ entry.row()  // J(i,2) =
                                 ) += entry.value() * x( entry.col() ); //        = c_21 * x_1


                // d(c_21 * x_2 * x_1) / dx_1 = c_21 * x_2
                jacobian.coeffRef( /* row == constr_id: */ constr_id
                                 , /* col == var_id:    */ entry.col()  // J(i,1) =
                                 ) += entry.value() * x( entry.row() ); //        = c_21 * x_2
            }

        } //...entries

    return jacobian;
} //...getJacobian

template <typename _Scalar> int
OptProblem<_Scalar>::precalcHessianCoeffs()
{
//    if ( _quadConstrList.size() )
//    {
//        std::cerr << "[" << __func__ << "]: " << "Hessian NOT implemented for quadratic constraints yet..." << std::endl;
//        throw new std::runtime_error( "[OptProblem::estimateHessianOfLagrangian] Hessian NOT implemented for quadratic constraints yet..." );
//    }

    _hessians.resize( /* obj */ 1 + /* constraints: */ this->getConstraintCount() );

    _hessians[0] = SparseMatrix( getVarCount(), getVarCount() );
    _hessians[0].reserve( _quadObjList.size() );
    for ( size_t i = 0; i != _quadObjList.size(); ++i )
    {
        // (c * x^2)'' == 2c. Second derivative has a 2 multiplier if it's a squared variable.
        if ( _quadObjList[i].row() != _quadObjList[i].col() )
            _hessians[0].coeffRef( _quadObjList[i].row(), _quadObjList[i].col() ) += _quadObjList[i].value();
        else
        {
            //std::cout << "[" << __func__ << "]: " << "inserting _hessians[0]" << _quadObjList[i].row() << ", " << _quadObjList[i].col() << " = " << _quadObjList[i].value() << std::endl;
            _hessians[0].coeffRef( _quadObjList[i].row(), _quadObjList[i].col() ) += _Scalar(2) * _quadObjList[i].value();
        }
    } // for linConstrList

    size_t max_j_nnz = 0;
    for ( size_t j = 0; j != this->getConstraintCount(); ++j )
    {
        _hessians[1+j] = SparseMatrix( getVarCount(), getVarCount() );

        if ( j < _quadConstrList.size() ) //if quadratic constraint exists
        {
            for ( size_t c = 0; c != _quadConstrList[j].size(); ++c )
            {
                SparseEntry const& entry = _quadConstrList[j][c];
                _hessians[1+j].coeffRef( entry.row(), entry.col() ) += entry.value();
                _hessians[1+j].coeffRef( entry.col(), entry.row() ) += entry.value();
            } //...for constraint entries
        } //...if quadratic constraint exists

        if ( _hessians[1+j].nonZeros() )
            max_j_nnz = 1+j;
    } //...for constraints

    if ( max_j_nnz < _hessians.size()-1 )
    {
        std::cout << "resizing hessians to " << max_j_nnz+1 << std::endl;
        _hessians.resize( max_j_nnz + 1 );
    }

    return EXIT_SUCCESS;
} // ...OptProblem::precalcHessianCoeffs()

template <typename _Scalar> typename OptProblem<_Scalar>::SparseMatrix const
OptProblem<_Scalar>::getHessian( Scalar const obj_factor, VectorX const& lambdas ) const
{
//    if ( lambdas.size() != _hessians.size() - 1 )
//    {
//        std::cerr <<"[" << __func__ << "]: " << "lambdas.size() " << lambdas.size() << " != " << _hessians.size() - 1 << " _hessians.size()-1" << std::endl;
//        throw new OptProblemException( "getHessian size error" );
//    }

    if ( _hessians.size() == 1 )
        return _hessians[0] * obj_factor;
    else
    {
        SparseMatrix hessian = _hessians[0] * obj_factor;
        for ( size_t j = 0; j != _hessians.size()-1; ++j )
        {
            hessian += _hessians[1+j] * lambdas[j];
        }
        return hessian;
    }
}

//______________________________________________________________________________

template <typename _Scalar> int
OptProblem<_Scalar>::addLinConstraint( BOUND bound_type, Scalar lower_bound, Scalar upper_bound, std::vector<Scalar> coeffs, LINEARITY c_lin  )
{
    // usage:
    //    coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n >= lower_bound              , bound_type == BOUND::GREATER_EQ
    // OR coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n <= upper_bound              , bound_type == BOUND::LESS_EQ
    // OR coeffs[0] * x_0 + coeffs[1] * x_1 ... + coeffs[n] * x_n =  lower_bound = upper_bound, bound_type == BOUND::FIXED
    if ( coeffs.size() != getVarCount() )
    {
        std::cerr << "[" << __func__ << "]: " << "A line in the constraints matrix A has to be varCount " << getVarCount() << " long, not " << coeffs.size() << std::endl;
        return EXIT_FAILURE;
    }

    // prepare
    SparseMatrix row_vector( 1, coeffs.size() );
    // add coeffs from new line
    for ( size_t col = 0; col != coeffs.size(); ++col )
    {
        // add non-zero elements to sparse representation
        if ( coeffs[col] != Scalar(0) )
        {
            row_vector.insert( 0, col ) = coeffs[col];
        }
    }

    // work
    return this->addConstraint( bound_type, lower_bound, upper_bound, &row_vector, c_lin );
}

template <typename _Scalar> int
OptProblem<_Scalar>::addConstraint( Variable const& constr )
{
    return this->addConstraint( /*   bound type: */ constr._bkx
                              , /*  lower bound: */ constr._blx
                              , /*  upper bound: */ constr._bux
                              , /* coeff_vector: */ NULL
                              , /*    linearity: */ constr._lin_x );
}

template <typename _Scalar> int
OptProblem<_Scalar>::addConstraint( BOUND bound_type, Scalar lower_bound, Scalar upper_bound, SparseMatrix *row_vector /* = NULL */, LINEARITY linearity /* = LINEAR */ )
{
    typedef typename SparseMatrix::Index Index;
    // check bounds
    if      ( bound_type == BOUND::EQUAL )
    {
        if ( lower_bound != upper_bound )
        {
            std::cerr << "[" << __func__ << "]: " << "[Warning] If bound_type is FIXED, upper_bound " << upper_bound << " == " << lower_bound << " lower_bound must hold" << std::endl;
            //return EXIT_FAILURE;
        }
    }

    // contraint matrix row
    const int row = _bkc.size();
    // var bound type
    _bkc.push_back( bound_type  );
    // lower bounds
    _blc.push_back( lower_bound );
    // upper bounds
    _buc.push_back( upper_bound );
    // constriant linearity
    _lin_c.push_back( linearity );

    // constraint coeffs
    if ( row_vector )
    {
        if ( (row_vector->cols() != static_cast<Index>(this->getVarCount())) || (row_vector->rows() != static_cast<Index>(1)) )
        {
            std::cerr << "[" << __func__ << "]: " << "row_vector->cols()(" << row_vector->cols() << ") != this->getVarCount() (" << this->getVarCount() << ") || row_vector->rows()( " << row_vector->rows() << " != 1, returning" << std::endl;
            return EXIT_FAILURE;
        }

        // add coeffs from new line
        for ( int k = 0; k != row_vector->outerSize(); ++k )
        {
            for ( typename SparseMatrix::InnerIterator it(*row_vector,k); it; ++it )
            {
                _linConstrList.push_back( SparseEntry(row, it.col(), it.value()) );
            }
        }
    }

    return EXIT_SUCCESS;
} // ...OptProblem::addConstraint

template <typename _Scalar> int
OptProblem<_Scalar>::addLinConstraints( SparseMatrix const& mx )
{
    typedef typename SparseMatrix::Index Index;

    if ( mx.rows() != static_cast<Index>(this->getConstraintCount()) || mx.cols() != static_cast<Index>(this->getVarCount()) )
    {
        std::cerr << "[" << __func__ << "]: " << "mx has to be m x n, m == constrCount(" << this->getConstraintCount() << "), n == varCount(" << this->getVarCount() << ")...returning" << std::endl;
        return EXIT_FAILURE;
    }

    this->_linConstrList.reserve( this->_linConstrList.size() + mx.nonZeros() );
    for ( int row = 0; row != mx.outerSize(); ++row )
        for ( typename SparseMatrix::InnerIterator it(mx,row); it; ++it )
        {
            this->_linConstrList.push_back( SparseEntry(it.row(), it.col(), it.value()) );
        } // ... for cols

    return EXIT_SUCCESS;
} // ...OptProblem::setLinConstraints()

template <typename _Scalar> int
OptProblem<_Scalar>::setLinConstraints( SparseMatrix const& mx )
{
    _linConstrList.clear();
    return this->addLinConstraints( mx );;
} // ...OptProblem::setLinConstraints()

template <typename _Scalar> int
OptProblem<_Scalar>::addQConstraint( size_t constr_id, size_t i, size_t j, Scalar coeff )
{
    if ( constr_id >= getConstraintCount() )
    {
        std::cerr << "constr_id <= getConstraintCount, please use addLinConstraint to initialize constr_id first! exiting" << std::endl;
        return EXIT_FAILURE;
    }

    if ( _quadConstrList.size() <= constr_id )
        _quadConstrList.resize( constr_id + 1 );

    // only lower triangle
    if ( j > i )
    {
        //std::cerr << "[" << __func__ << "]: " << "swapping i and j to make it lower triangular" << std::endl;
        std::swap( i, j );
    }

    _quadConstrList[constr_id].push_back( SparseEntry(i,j,coeff) );
    //if ( verbose ) std::cout << "[" << __func__ << "]: " << "qconstrlist is now " << _quadConstrList.size() << "( " << _quadConstrList[0].size() << ")" << " long" << std::endl;

    return EXIT_SUCCESS;
} // ...OptProblem::addQConstraint

template <typename _Scalar> int
OptProblem<_Scalar>::addQConstraints( SparseMatrix const& mx )
{
    typedef typename SparseMatrix::Index Index;

    if ( mx.rows() != static_cast<Index>(this->getVarCount()) || mx.cols() != static_cast<Index>(this->getVarCount()) )
    {
        std::cerr << "[" << __func__ << "]: " << "mx has to be n x n, where n == varCount(" << this->getVarCount() << ")...returning" << std::endl;
        return EXIT_FAILURE;
    }
    if ( this->_quadConstrList.size() >= this->getConstraintCount() )
    {
        std::cerr << "[" << __func__ << "]: " << "this->_quadConstrList.size() " << this->_quadConstrList.size() << " >= " << this->getConstraintCount() << " this->getConstraintCount(), call addLinConstr first!" << std::endl;
        throw new OptProblemException( "Too many quadratic constraints" );
    }

    this->_quadConstrList.push_back( SparseEntries() );
    this->_quadConstrList.back().reserve( mx.nonZeros() );
    for ( int row = 0; row != mx.outerSize(); ++row )
        for ( typename SparseMatrix::InnerIterator it(mx,row); it; ++it )
        {
            this->_quadConstrList.back().push_back( SparseEntry(it.row(), it.col(), it.value()) );
        } // ... for cols

    return EXIT_SUCCESS;
} // ...OptProblem::setLinConstraints()

template <typename _Scalar> typename OptProblem<_Scalar>::SparseMatrix
OptProblem<_Scalar>::getLinConstraintsMatrix() const
{
    SparseMatrix mx( this->getConstraintCount(), this->getVarCount() );
    mx.setFromTriplets( this->_linConstrList.begin(), this->_linConstrList.end() );

    return mx;
} // ...OptProblem::getLinConstraintsMatrix()

template <typename _Scalar> typename OptProblem<_Scalar>::SparseMatrix
OptProblem<_Scalar>::getQuadraticConstraintsMatrix( size_t j ) const
{
    SparseMatrix mx( this->getVarCount(), this->getVarCount() );

    if ( this->getQuadraticConstraints().size() > j )
    {
        mx.setFromTriplets( this->getQuadraticConstraints(j).begin(), this->getQuadraticConstraints(j).end() );
    }

    return mx;
} // ...OptProblem::getLinConstraintsMatrix()

template <typename _Scalar> int
OptProblem<_Scalar>::setStartingPoint( typename OptProblem<_Scalar>::SparseMatrix const& x0 )
{
    int err = EXIT_SUCCESS;

    if ( (x0.outerSize() != 1) && (x0.innerSize() != 1) )
    {
        std::cerr << "[" << __func__ << "]: " << "x0 has to be a vector! It's now " << x0.outerSize() << "x" << x0.innerSize() << std::endl;
        err = EXIT_FAILURE;
    }

    if ( EXIT_SUCCESS == err )
    {
        _x0.resize( x0.rows(), x0.cols() );
        _x0.setZero();

        for ( int k = 0; k < x0.outerSize(); ++k )
            for ( typename OptProblem<_Scalar>::SparseMatrix::InnerIterator it(x0, k); it; ++it )
            {
                _x0( it.row(), it.col() ) = it.value();
            }

        _useStartingPoint = true;
    }

    return err;
} //...OptProblem::setStartingPoint()

template <typename _Scalar> typename OptProblem<_Scalar>::SparseMatrix
OptProblem<_Scalar>::getStartingPointMatrix() const
{
    SparseMatrix mx( this->_x0.rows(), 1 );
    for ( int r = 0; r != this->_x0.rows(); ++r )
        mx.insert( r, 0 ) = this->_x0( r );
    return mx;
} // ...OptProblem::getLinConstraintsMatrix()

//______________________________________________________________________________

template <typename _SparseMatrixT>
int _printSparseMatrix( _SparseMatrixT const& mx, std::string const& title = "smx", const int entry_limit = 10 )
{
    int cnt = 0;
    std::cout << title << "(" << mx.nonZeros() << " entries): ";
    if ( mx.nonZeros() )
    {
        for ( int row = 0; (row != mx.outerSize()) && (cnt < entry_limit); ++row )
        {
            for ( typename _SparseMatrixT::InnerIterator it(mx,row); it && (cnt < entry_limit); ++it, ++cnt )
            {
                std::cout << "(" << it.row() << ", " << it.col() << "," << it.value() << "), ";
            }
        }
        if ( cnt < mx.nonZeros() )
            std::cout << "..."; //and " << mx.nonZeros() - cnt << " other entries...";
    }
    else
        std::cout << "empty";

    std::cout << std::endl;
    return EXIT_SUCCESS;
}

template <typename _Scalar> int
OptProblem<_Scalar>::printProblem( size_t entry_limit /* = 10 */ ) const
{
    // Const objective
    std::cout << "[" << __func__ << "]: " << "const objective: " << this->getObjectiveBias() << std::endl;

    // Linear objectives
    {
        std::cout<< "[" << __func__ << "]: " << "_linObjs(" <<_linObjs.size()<<" entries): ";
        for ( size_t vi=0; (vi!=_linObjs.size()) && (vi < entry_limit); ++vi )
            std::cout << _linObjs[vi] <<" ";
        if ( _linObjs.size() > entry_limit )
            std::cout << "...";
        std::cout << "\n";
    }

    // Quad obj
    {
        SparseMatrix Qo = this->getQuadraticObjectivesMatrix();
        std::cout << "[" << __func__ << "]: ";
        _printSparseMatrix( Qo, "Qo", entry_limit );
    }

    // Linear constraints
    {
        SparseMatrix A = this->getLinConstraintsMatrix();
        std::cout << "[" << __func__ << "]: ";
        _printSparseMatrix( A, "A", entry_limit );
    }

    // Quadratic constraints
    if ( this->getConstraintCount() <= entry_limit )
    {
        for ( size_t i = 0; i != this->getConstraintCount(); ++i )
        {
            SparseMatrix Qi = this->getQuadraticConstraintsMatrix( i );
            if ( Qi.nonZeros() )
            {
                char name[255]; sprintf(name, "Q%lu", i);
                _printSparseMatrix( Qi, name, entry_limit );
            }
        }
    }
    else
    {
        size_t cnt = 0;
        for ( size_t i = 0; i != this->getConstraintCount(); ++i )
        {
            cnt += ( !this->getQuadraticConstraints(i).empty() );
        }
        std::cout << "[" << __func__ << "]: " << " have " << cnt << " nonempty quadratic constraint matrices" << std::endl;
    }


    return EXIT_SUCCESS;
} // ...OptProblem::printProblem()

//_____________________________________________________________________________

template <typename _Scalar> int
OptProblem<_Scalar>::write( std::string const& path ) const
{
    const int entry_limit = 5;
    int entries = 0;

    std::cout << "[" << __func__ << "]: " << "creating " << path << std::endl;
    mkdir( path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    std::vector<std::string> paths;

    // vars
    paths.push_back( getAuxName() );
    {
        std::ofstream aux_file( path + "/" + paths.back() );
        if ( !aux_file.is_open() )
        {
            std::cerr << "[" << __func__ << "]: " << "could not open " << path + "/" + paths.back() << std::endl;
            return EXIT_FAILURE;
        }

        aux_file << "# vars, constraints, objective bias (c)" << std::endl;
        aux_file << this->getVarCount() << "," << this->getConstraintCount() << "," << this->getObjectiveBias() << std::endl;
        aux_file << "# bound_type, var_type, lower_boud, upper_bound, linearity" << std::endl;
        for ( size_t j = 0; j != this->getVarCount(); ++j )
        {
            aux_file << this->getVarBoundType(j) << ","
                     << this->getVarType( j ) << ","
                     << this->getVarLowerBound( j ) << ","
                     << this->getVarUpperBound( j ) << ","
                     << this->getVarLinearity( j )
                     << std::endl;
        } // for each variable
        aux_file << "# bound_type, lower_boud, upper_bound, linearity" << std::endl;
        for ( size_t i = 0; i != this->getConstraintCount(); ++i )
        {
            aux_file << this->getConstraintBoundType( i ) << ","
                     << this->getConstraintLowerBound( i ) << ","
                     << this->getConstraintUpperBound( i ) << ","
                     << this->getConstraintLinearity( i )
                     << std::endl;
        } // for each variable

        std::cout << "[" << __func__ << "]: " << "wrote " << path + "/" + paths.back() << std::endl;
    } // vars file

    // qo
    paths.push_back(getqoName());
    io::writeSparseMatrix( this->getLinObjectivesMatrix(), path + "/" + paths.back(), 0 );
    std::cout << "[" << __func__ << "]: " << "wrote " << path + "/" + paths.back() << "\t...";

    // Qo
    paths.push_back(getQoName());
    io::writeSparseMatrix( this->getQuadraticObjectivesMatrix(), path + "/" + paths.back(), 0 );
    std::cout << "[" << __func__ << "]: " << "wrote " << path + "/" + paths.back() << "\t...";

    // A
    paths.push_back(getAName());
    io::writeSparseMatrix( this->getLinConstraintsMatrix(), path + "/" + paths.back(), 0 );
    std::cout << "[" << __func__ << "]: " << "wrote " << path + "/" + paths.back() << "\t...";

    // Qi
    for ( size_t i = 0; i != this->getQuadraticConstraints().size(); ++i )
    {
        paths.push_back( getQiName(i) );
        io::writeSparseMatrix( this->getQuadraticConstraintsMatrix(i), path + "/" + paths.back(), 0 );
        if ( entries++ < entry_limit ) std::cout << "[" << __func__ << "]: " << "wrote " << path + "/" + paths.back() << "\t...";
    }

    // X0
    paths.push_back( getX0Name() );
    io::writeSparseMatrix( this->getStartingPointMatrix(), path + "/" + paths.back(), 0 );
    std::cout << "[" << __func__ << "]: " << "wrote " << path + "/" + paths.back() << "\t...";

    std::string proj_path = path + "/problem.proj";
    std::ofstream proj_file( proj_path );
    if ( !proj_file.is_open() )
    {
        std::cerr << "[" << __func__ << "]: " << "could not open " << proj_path << std::endl;
    }
    for ( size_t i = 0; i != paths.size(); ++i )
        proj_file << paths[i] << std::endl;

    proj_file.close();
    std::cout << "[" << __func__ << "]: " << "project written to " << proj_path << std::endl;

    return EXIT_FAILURE;
} // ... Optproblem::write

template <typename _Scalar> int
OptProblem<_Scalar>::_parseAuxFile( std::string const& aux_file_path )
{
    int err = EXIT_SUCCESS;

    // open file
    std::ifstream aux_file( aux_file_path );
    if ( !aux_file.is_open() )
    {
        std::cerr << "[" << __func__ << "]: " << "Coult not open aux_file " << aux_file_path << std::endl;
        return EXIT_FAILURE;
    }

    // parse lines
    int         vars        = -1, // "uninited"
                constraints = -1; // "uninited"
    int         lid         = -1; // -1 == there's an extra line to read first, the header
    std::string line;             // tmp storage of read line
    while ( getline(aux_file, line) )
    {
        // skip comment
        if ( line[0] == '#') continue;

        // parse line
        std::istringstream iss( line );
        std::string        word;
        int                word_id = 0; // line parse state
        Variable           tmp;         // output storage variable
        while ( std::getline(iss, word, ',') )
        {
            if ( vars < 0 ) // header line
            {
                vars = atoi( word.c_str() );        // parse number of variables
                std::getline(iss, word, ',');
                constraints = atoi( word.c_str() ); // parse number of constraints
                std::getline(iss, word, ',');
                this->_cfix = atof( word.c_str() ); // parse constant objective function bias

                std::cout << "bias is now " << this->getObjectiveBias()
                          << " reading " << vars << " vars and " << constraints << " constraints"
                          << std::endl;

                std::getline( iss, word );          // clear rest of line
            }
            else // variable or constraint line
            {
                if ( lid < vars ) // variable line
                {
                    switch ( word_id ) // bound_type, var_type, lower_boud, upper_bound, linearity
                    {
                        case 0: // Variable bound type
                            tmp._bkx = static_cast<BOUND>( atoi(word.c_str()) );
                            break;
                        case 1: // Variable type (cont/integer/binary)
                            tmp._type_x = static_cast<VAR_TYPE>( atoi(word.c_str()) );
                            break;
                        case 2: // Variable lower bound
                            tmp._blx = atof( word.c_str() );
                            break;
                        case 3: // Variable upper bound
                            tmp._bux = atof( word.c_str() );
                            break;
                        case 4: // Variable linearity
                            tmp._lin_x = static_cast<LINEARITY>( atoi(word.c_str()) );
                            break;
                        default:
                            std::cerr << "[" << __func__ << "]: " << "Unknown VAR word_id..." << word_id << std::endl;
                            break;
                    } //...switch word_id
                } //...if variable
                else // constraint line
                {
                    switch( word_id ) // bound_type, lower_boud, upper_bound, linearity
                    {
                        case 0: // Constraint bound type
                            tmp._bkx = static_cast<BOUND>( atoi(word.c_str()) );
                            break;
                        case 1: // Constraint lower bound
                            tmp._blx = atof( word.c_str() );
                            break;
                        case 2: // Constraint upper bound
                            tmp._bux = atof( word.c_str() );
                            break;
                        case 3: // Constraint linearity
                            tmp._lin_x = static_cast<LINEARITY>( atoi(word.c_str()) );
                            break;
                        default:
                            std::cerr << "[" << __func__ << "]: " << "Unknown CNSTR word_id..." << word_id << ", word: " << word << std::endl;
                            break;
                    } //...switch word_id
                } //...if constraint

                // increment line parse status
                ++word_id;
            } //... non-header line
        } //...while line has more words

        // add variable/constraint
        if ( lid >= 0 )
        {
            if ( lid < vars ) this->addVariable  ( tmp ); // variable
            else              this->addConstraint( tmp ); // constraint
        } // if non-header-line

        // increase read line count (did not count comments)
        ++lid;
    } //... while file has more lines

    // close file
    aux_file.close();

    return err;
} //...OptProblem::_parseAuxFile()

template <typename _Scalar> int
OptProblem<_Scalar>::read( std::string proj_file_path )
{
    int err = EXIT_SUCCESS;
    const int entry_limit = 5;
    int entries = 0;

    // Parse project path
    std::string proj_path = ".";
    if ( proj_file_path.rfind(".proj") != proj_file_path.size() - 5 )
    {
        proj_file_path = proj_file_path + "/problem.proj";
    }

    size_t slash_pos = proj_file_path.rfind("/");
    if ( slash_pos != std::string::npos )
    {
        proj_path = proj_file_path.substr( 0, slash_pos );
    }

    // Open project file
    std::vector<std::string> paths;
    {
        // open
        std::ifstream proj_file( proj_file_path );
        if ( !proj_file.is_open() )
        {
            std::cerr << "[" << __func__ << "]: " << "could not open " << proj_file_path << std::endl;
            return EXIT_FAILURE;
        }

        // parse project file
        std::string line;
        while ( getline(proj_file, line) )
        {
            if ( line[0] == '#') continue;

            paths.push_back( line );
        }
        proj_file.close();
    } //... parse project file


    // Parse project files
    for ( size_t i = 0; i != paths.size(); ++i )
    {
        std::string fname = paths[i];
        if ( entries < entry_limit ) std::cout << "[" << __func__ << "]: " << "reading " << proj_path + "/" + paths[i] << "..."; fflush(stdout);
        // aux file or sparse matrix
        if ( fname.find("aux") != std::string::npos )
        {
            this->_parseAuxFile( proj_path + "/" + paths[i] );
        }
        else // sparse matrix
        {
            // read
            SparseMatrix mx = io::readSparseMatrix<_Scalar>( proj_path + "/" + paths[i], 0 );

            // save
            if ( fname.find("qo") != std::string::npos ) // lin objective
            {
                std::cout << "read " << proj_path + "/" + paths[i] << " as lin obj mx" << ",  ";
                err += this->addLinObjectives( mx );
                std::cout << "problem varcount : " << this->getVarCount() << std::endl;
            } //...if qo
            else if ( fname.find("Qo") != std::string::npos ) // quadratic objective
            {
                std::cout << "read " << proj_path + "/" + paths[i] << " as quadratic obj mx" << ", ";
                err += this->addQObjectives( mx );
            } //...if Qo
            else if ( fname.find("A") != std::string::npos ) // lin constraint
            {
                std::cout << "read " << proj_path + "/" + paths[i] << " as linear constraints mx" << ", ";
                err += this->addLinConstraints( mx );
            } //...if A
            else if ( fname.find("Q") != std::string::npos ) // Qo will be parsed before, so that's fine
            {
                if ( entries++ < entry_limit ) std::cout << "read " << proj_path + "/" + paths[i] << " as quadratic constraints mx" << ", ";
                err += this->addQConstraints( mx );
            }
            else if ( fname.find(this->getX0Name()) != std::string::npos ) // Starting point
            {
                std::cout << "read " << proj_path + "/" + paths[i] << " as starting point mx" << ", ";
                if ( mx.nonZeros() )
                    err += this->setStartingPoint( mx );
                else
                    std::cerr << "[" << __func__ << "]: " << "X0 has no entries, not setting starting point" << std::endl;
            } //...ifelse sparse matrix type
        } //...sparse matrix
    } //...for project file
    if ( entries < entry_limit ) std::cout << std::endl;

    // print summary
    this->printProblem();

    return err;
} // ...OptProblem::read

} // ...namespace qcqpp

#endif // QCQPCPP_SGOPTPROBLEM_HPP
