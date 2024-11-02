#define UL_TRIANG_MATRIX_CPP
#include "ul_triangMatrix.hpp"

template <typename T>
const ul_triangMatrix<T> ul_triangMatrix<T>::staticHelper;

template <typename T>
T &ul_triangMatrix<T>::operator()(const size_t row, const size_t col)
{
    if (row > col)
        return this->_begin[(((row - 1) * row) >> 1) + col];
    throw "index inaccessable";
    return this->_begin[-1];
};

template <typename T>
template <typename U>
ul_triangMatrix<T> *ul_triangMatrix<T>::holdInv(const ul_triangMatrix<U> &other, const bool checkSize)
{
    if (checkSize)
        MatrixBase<T>::resizeLike(other, false, false);

    for (size_t i = 0; i < this->_rows; i++)
    {
        const size_t other_index = (((i - 1) * i) >> 1);
        for (size_t j = 0; j < i; j++)
        {
            T sum = other._begin[other_index + j];

            for (size_t k = j + 1; k < i; k++)
                sum += other._begin[other_index + k] * this->_begin[(((k - 1) * k) >> 1) + j];

            this->_begin[other_index + j] = -sum;
        }
    }

    return this;
};

// ul_triangMatrix and ul_triangMatrix
template <typename T>
template <typename U, typename V>
ul_triangMatrix<T> *ul_triangMatrix<T>::holdMul(const ul_triangMatrix<U> &a, const ul_triangMatrix<V> &b, const bool checkSize, const bool checkOverlap)
{
    if (checkSize)
        MatrixBase<T>::checkSize_mul(a, b);
    if (checkOverlap)
        MatrixBase<T>::checkOverlap(a, b);

    size_t index = 0;
    for (size_t i = 0; i < this->_rows; i++)
    {
        const size_t a_index = (((i - 1) * i) >> 1);
        for (size_t j = 0; j < i; j++)
        {
            T sum = a._begin[a_index + j] + b._begin[a_index + j];

            for (size_t k = j + 1; k < i; k++)
                sum += a._begin[a_index + k] * b._begin[(((k - 1) * k) >> 1) + j];

            this->_begin[index++] = sum;
        }
    }
    return this;
};