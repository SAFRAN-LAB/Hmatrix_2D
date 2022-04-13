#ifndef __HODLR_Matrix__
#define __HODLR_Matrix__

#ifdef MKL_ENABLED
    #define EIGEN_USE_MKL_ALL
#endif

#include "myHeaders.hpp"

class HODLR_Matrix
{
public:

// Size of the matrix:
int N;

// Modulo operator:
// This is separately defined to make sure
// that positive values are always returned
int mod(int a, int b)
{
        return ((a % b + b) % b);
}

// Constructor:
explicit HODLR_Matrix(int N)
{
        this->N = N;
}

// Returns individual entries of the matrix:
virtual dtype getMatrixEntry(int j, int k)
{
        // FROM EXPERIENCE: Incase the user makes a mistake in
        // setting the derived class, this warns the user:
        std::cout << "Returning zero! Ensure that derived class is properly set!" << std::endl;
        return 0.0;
}

Vec getRow(int j, int n_col_start, int n_cols);
Vec getCol(int k, int n_row_start, int n_rows);
Vec getDiag1(int j, int k, int n_rows, int n_cols);
Vec getDiag2(int j, int k, int n_rows, int n_cols);
Mat getMatrix(int j, int k, int n_rows, int n_cols);

/* / Destructor:
   ~HODLR_Matrix() {
   };*/
virtual ~HODLR_Matrix() {
}
};

#endif /*__HODLR_Matrix__*/
