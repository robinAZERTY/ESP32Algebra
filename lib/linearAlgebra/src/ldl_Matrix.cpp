#define LDL_MATRIX_CPP
#include "ldl_Matrix.hpp"

template <typename T>
ldl_matrix<T> *ldl_matrix<T>::decompose()
{
    for (size_t j = 0; j < this->_rows; j++)
    {
        D[j] = this->_begin[(j * (j + 1) >> 1) + j];
        const size_t l_index = ((j - 1) * j) >> 1;
        for (size_t k = 0; k < j; k++)
            D[j] -= L[l_index + k] * L[l_index + k] * D[k];
        for(size_t i = j + 1; i < this->_rows; i++)
        {
            const size_t l_index2 = ((i - 1) * i) >> 1;
            L[l_index2 + j] = this->_begin[(i * (i + 1) >> 1) + j];
            for (size_t k = 0; k < j; k++)
                L[l_index2 + j] -= L[l_index + k] * L[l_index2 + k] * D[k];
            L[l_index2 + j] /= D[j];
        }
    }
    return referT();
}