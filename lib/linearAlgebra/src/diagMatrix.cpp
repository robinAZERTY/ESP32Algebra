#define DIAG_MATRIX_CPP
#include "diagMatrix.hpp"

template <typename T>
const diagMatrix<T> diagMatrix<T>::staticHelper;

template <typename T>
T &diagMatrix<T>::operator()(const size_t row, const size_t col)
{
    if (row == col)
        return this->_begin[row];
    else
        throw "diagMatrix::operator() out of range";
}

template <typename T>
T diagMatrix<T>::det() const
{
    if (this->rows() != this->cols())
    {
        throw "diagMatrix::det() not a square matrix";
        return 0;
    }
    T res = 1;
    for (size_t i = 0; i < this->_cols; i++)
        res *= this->_begin[i];
    return res;
}