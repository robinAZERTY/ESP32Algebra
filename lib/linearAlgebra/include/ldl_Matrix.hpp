#ifndef LDL_MATRIX_HPP
#define LDL_MATRIX_HPP

#include "ul_triangMatrix.hpp"
#include "diagMatrix.hpp"
#include "uu_triangMatrix.hpp"
// #include "triangMatrix.hpp"

template <typename T>
class ldl_matrix : public symMatrix<T>
{
    // protected:
    public:
        ldl_matrix<T> *referT() noexcept {LT.refer(L, this->_cols, this->_rows); return this;}
    public:
        ul_triangMatrix<T> L;
        diagMatrix<T> D;
        uu_triangMatrix<T> LT;

        ldl_matrix() : symMatrix<T>(){}
        ldl_matrix(const size_t order) : symMatrix<T>(){this->resize(order);}
        ldl_matrix<T> *resize(const size_t order, const bool deallocIfPossible = false, const bool saveData = true) {symMatrix<T>::resize(order,order, deallocIfPossible, saveData); L.resize(order,order); D.resize(order,order); referT(); return this;}
        ldl_matrix<T> *decompose();
};


#ifndef LDL_MATRIX_CPP
#include "ldl_Matrix.cpp"
#endif
#endif