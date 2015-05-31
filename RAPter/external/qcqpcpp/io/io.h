#ifndef QCQPCPP_IO_H
#define QCQPCPP_IO_H

#include <fstream>
#include <iomanip>

#include "Eigen/Sparse"

namespace qcqpcpp {
namespace io {

//! \brief readSparseMatrix Reads sparse matrix text file. First row: rows,cols, then in each line: row,col,value triplets
//! \param path             Path of txt file containing the sparse matrix
//! \param index_offset     Adds this to every row and col value in the triplets. Set to -1 to convert from MATLAB dumps.
//! \param return           RowMajor Sparsematrix containing the values
template <typename Scalar> inline Eigen::SparseMatrix<Scalar,Eigen::RowMajor>
readSparseMatrix( std::string path, int index_offset /* = -1 */ )
{
    Eigen::SparseMatrix<Scalar,Eigen::RowMajor> mx;
    std::vector< Eigen::Triplet<Scalar> > entries;

    std::ifstream f_matrix( path );
    if ( !f_matrix.is_open() )
    {
        std::cerr << "[" << __func__ << "]: " << "could not open " << path << "...exiting\n";
        return mx;
    }

    int rows = -1, cols = -1;
    std::string line;
    while ( getline(f_matrix, line) )
    {
        // skip comments
        if ( line[0] == '#') continue;

        // read word by word
        std::istringstream  iss( line );
        std::string         word;
        int                 word_id = 0;
        int i = -1, j = -1;
        Scalar val;
        while ( std::getline(iss, word, ',') )
        {
            if ( rows < 0 )
            {
                rows = atoi( word.c_str() );
                std::getline(iss, word, ',');
                cols = atoi( word.c_str() );
                std::getline(iss, word, ',');
            }
            else switch ( word_id )
            {
                case 0:
                    i = atoi( word.c_str() ) + index_offset;
                    //if ( i >= rows ) rows = i + 1;
                    break;
                case 1:
                    j = atoi( word.c_str() ) + index_offset;
                    //if ( j >= cols ) cols = j + 1;
                    break;
                case 2:
                    val = atof( word.c_str() );
                    break;
                default:
                    std::cout << "[" << __func__ << "]: " << "word_id not 0,1,2?" << std::endl;
                    break;
            }
            ++word_id;
        } // ... while getword

        if ( (i >= 0) && (j>=0) )
            entries.push_back( Eigen::Triplet<Scalar>(i,j,val) );
    } // ... while getline

    mx.resize( rows, cols );
    mx.setFromTriplets( entries.begin(), entries.end() );

    f_matrix.close();

    return mx;
} // ... readSparseMatrix

template <typename Scalar> inline int
writeSparseMatrix( Eigen::SparseMatrix<Scalar,Eigen::RowMajor> const& mx, std::string path, int index_offset )
{
    std::ofstream f_matrix( path );
    if ( !f_matrix.is_open() )
    {
        std::cerr << "[" << __func__ << "]: " << "could not open " << path << "...exiting\n";
        return EXIT_FAILURE;
    }

    f_matrix << mx.rows() << "," << mx.cols() << std::endl;
    for ( int row = 0; row != mx.outerSize(); ++row )
        for ( typename Eigen::SparseMatrix<Scalar,Eigen::RowMajor>::InnerIterator it(mx,row); it; ++it )
            f_matrix << it.row() + index_offset << ","
                     << it.col() + index_offset << ","
                     << std::setprecision(16) << it.value()
                     << std::endl;

    f_matrix.close();
    return EXIT_SUCCESS;
} //... writeSparseMatrix

} //...namespace io
} //...namespace qcqpcpp

#endif // QCQPCPP_IO_H
