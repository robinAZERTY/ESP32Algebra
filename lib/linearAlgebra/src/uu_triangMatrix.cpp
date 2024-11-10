#define UU_TRIANG_MATRIX_CPP
#include "uu_triangMatrix.hpp"

template <typename T>
const uu_triangMatrix<T> uu_triangMatrix<T>::staticHelper;

template <typename T>
T &uu_triangMatrix<T>::operator()(const size_t row, const size_t col)
{
    if (col > row)
        return this->_begin[(((col - 1) * col) >> 1) + row];
    throw "index inaccessable";
    return this->_begin[-1];
};
