#define UKF_CPP
#include "ukf.hpp"

#define max(a, b) (a > b ? a : b)

template <size_t x_dim, size_t u_dim, size_t c_dim, size_t z_num, typename T>
void Ukf<x_dim, u_dim, c_dim, z_num, T>::compute_S(const size_t z_idx)
{
}

template <size_t x_dim, size_t u_dim, size_t c_dim, size_t z_num, typename T>

Ukf<x_dim, u_dim, c_dim, z_num, T>::Ukf()
{
    X_sp.resize(2 * x_dim + 2 * u_dim, x_dim);
    X_sp.fill(0);
}

template <size_t x_dim, size_t u_dim, size_t c_dim, size_t z_num, typename T>

void Ukf<x_dim, u_dim, c_dim, z_num, T>::setMeasurementFunction(Vector_f2<T> h, size_t z_dim, size_t z_idx)
{
    this->h[z_idx] = h;
    this->z_dim[z_idx] = z_dim;
    this->S[z_idx].resize(z_dim);
    this->h_sp[z_idx].resize(2 * x_dim , z_dim);
    this->h_sp[z_idx].fill(0);
    this->h_val[z_idx].resize(z_dim);
    this->h_val[z_idx].fill(0);
}

template <size_t x_dim, size_t u_dim, size_t c_dim, size_t z_num, typename T>

void Ukf<x_dim, u_dim, c_dim, z_num, T>::compute_state_sigma_points()
{
    if (f == nullptr)
        throw "Ukf::predict() state transition function not set";
    this->P.decompose();
    this->Q.decompose();
    // compute state sigma points
    Vector<T> sp; // temporary sigma point
    for (size_t i = 0; i < x_dim; i++)
        this->P.D[i] = sqrt(this->P.D[i])*sqrt_x_dim_u_dim;
    for (size_t i = 0; i < u_dim; i++)
        this->Q.D[i] = sqrt(this->Q.D[i])*sqrt_x_dim_u_dim;

    internal::tmp<triangMatrix<T>> *tmp = internal::tmp<triangMatrix<T>>::get(max(x_dim, u_dim));
    // internal::tmp<Vector<T>> *sp = internal::tmp<Vector<T>>::get(x_dim); // temporary sigma point
    tmp->holdMul(this->P.L, this->P.D);
    for (size_t i = 0; i < x_dim; i++)
    {
        for (size_t j = i; j < x_dim; j++)
            this->X[j] += (*tmp)(j, i);
        sp = f(this->X, this->U, this->C);
        for (size_t j = 0; j < x_dim; j++)
            X_sp(i, j) = sp[j];
        for (size_t j = i; j < x_dim; j++)
            this->X[j] -= 2*(*tmp)(j, i);
        sp = f(this->X, this->U, this->C);
        for (size_t j = 0; j < x_dim; j++)
            X_sp(i + x_dim, j) = sp[j];
        for (size_t j = i; j < x_dim; j++)
            this->X[j] += (*tmp)(j, i);      
    }
    tmp->holdMul(this->Q.L, this->Q.D);
    for (size_t i = 0; i < u_dim; i++)
    {
        for (size_t j = i; j < u_dim; j++)
            this->U[j] += (*tmp)(j, i);
        sp = f(this->X, this->U, this->C);
        for (size_t j = 0; j < x_dim; j++)
            X_sp(i + 2*x_dim, j) = sp[j];
        for (size_t j = i; j < u_dim; j++)
            this->U[j] -= 2*(*tmp)(j, i);
        sp = f(this->X, this->U, this->C);
        for (size_t j = 0; j < x_dim; j++)
            X_sp(i + 2*x_dim + u_dim, j) = sp[j];
        for (size_t j = i; j < u_dim; j++)
            this->U[j] += (*tmp)(j, i);      
    }    
    tmp->release();
    // sp->release();
}

template <size_t x_dim, size_t u_dim, size_t c_dim, size_t z_num, typename T>

void Ukf<x_dim, u_dim, c_dim, z_num, T>::compute_measurement_sigma_points(const size_t z_idx)
{
    if (h[z_idx] == nullptr)
        throw "Ukf::update() measurement function not set";
    this->P.decompose();
    // compute state sigma points
    // internal::tmp<Vector<T>> *sp = internal::tmp<Vector<T>>::get(z_dim[z_idx]); // temporary sigma point
    Vector<T> sp; // temporary sigma point
    for (size_t i = 0; i < x_dim; i++)
        this->P.D[i] = sqrt(this->P.D[i])*sqrt_x_dim_u_dim;
    internal::tmp<triangMatrix<T>> *tmp = internal::tmp<triangMatrix<T>>::get(x_dim);
    tmp->holdMul(this->P.L, this->P.D);
    h_sp[z_idx].fill(-1);
    // compute measurement sigma points
    for (size_t i = 0; i < x_dim; i++)
    {
        for (size_t j = i; j < x_dim; j++)
            this->X[j] += (*tmp)(j, i);
            // this->X[j] += tmp->_begin[((j * (j + 1)) >> 1) + i];
        sp = h[z_idx](this->X, this->C);
        for (size_t j = 0; j < z_dim[z_idx]; j++)
            h_sp[z_idx](i, j) = sp[j];
        for (size_t j = i; j < x_dim; j++)
            this->X[j] -= 2*(*tmp)(j, i);
        sp = h[z_idx](this->X, this->C);
        for (size_t j = 0; j < z_dim[z_idx]; j++)
            h_sp[z_idx](i + x_dim, j) = sp[j];
        for (size_t j = i; j < x_dim; j++)
            this->X[j] += (*tmp)(j, i);
    }
    // sp->release();
    tmp->release();
}

template <size_t x_dim, size_t u_dim, size_t c_dim, size_t z_num, typename T>

void Ukf<x_dim, u_dim, c_dim, z_num, T>::predict()
{
    compute_state_sigma_points();
    // compute predicted state (mean of sigma points)
    X.fill(0);
    for (size_t i = 0; i < X_sp.rows(); i++)
        for (size_t j = 0; j < x_dim; j++)
            X[j] += X_sp(i, j);
    X /= X_sp.rows();
    // compute predicted state covariance
    P.fill(0);
    // for (size_t i = 0; i < X_sp.rows(); i++)
    //     for (size_t j = 0; j < x_dim; j++)
    //         for (size_t k = 0; k <=j; k++)
    //             P(j, k) += (X_sp(i, j) - X[j]) * (X_sp(i, k) - X[k]);

    for (size_t i = 0; i <  X_sp.rows(); i++)
        for (size_t k = 0; k < x_dim; k++)
            for (size_t j = k; j < x_dim; j++)
                P(j, k) += (X_sp(i, j) - X[j]) * (X_sp(i, k) - X[k]);
    P /= X_sp.rows();
}

template <size_t x_dim, size_t u_dim, size_t c_dim, size_t z_num, typename T>

void Ukf<x_dim, u_dim, c_dim, z_num, T>::update(const Vector<T> &Z, const symMatrix<T> &R, const size_t z_idx)
{
    if (z_dim[z_idx] != Z.size())
        throw "Ukf::update() measurement dimension mismatch";
    // compute_state_sigma_points();
    compute_measurement_sigma_points(z_idx);
    // compute predicted measurement (mean of sigma points)
    h_val[z_idx].fill(0);
    for (size_t i = 0; i < h_sp[z_idx].rows(); i++)
        for (size_t j = 0; j < z_dim[z_idx]; j++)
            h_val[z_idx][j] += h_sp[z_idx](i, j);
    h_val[z_idx] /= h_sp[z_idx].rows();
    // compute predicted measurement covariance
    S[z_idx].fill(0);
    for (size_t i = 0; i < h_sp[z_idx].rows(); i++)
        for (size_t j = 0; j < z_dim[z_idx]; j++)
            for (size_t k = 0; k <=j; k++)
                S[z_idx](j, k) += (h_sp[z_idx](i, j) - h_val[z_idx][j]) * (h_sp[z_idx](i, k) - h_val[z_idx][k]);
    S[z_idx] /= h_sp[z_idx].rows();
    S[z_idx] += R;
    S[z_idx].holdInv(S[z_idx], false);
    // compute cross covariance
    internal::tmp<rowMajorMatrix<T>> *Pxz = internal::tmp<rowMajorMatrix<T>>::get(x_dim, z_dim[z_idx]);
    Pxz->fill(0);
      for (size_t i = 0; i < 2*x_dim; i++)
        for (size_t j = 0; j < x_dim; j++)
            for (size_t k = 0; k <z_dim[z_idx]; k++)
                (*Pxz)(j, k) += (X_sp(i, j) - X[j]) * (h_sp[z_idx](i, k) - h_val[z_idx][k]);
    
    *Pxz /= 2*(int)x_dim;
    // K.resize(x_dim, z_dim[z_idx],false,true);
    K.holdMul((*Pxz), S[z_idx],true);
    K.referT();
    // update state
    // X.addMul(K,*((h_val[z_idx]-Z).release()));
    X.addMul(K,*((Z-h_val[z_idx]).release()));
    // update state covariance
    // P.subMul(*Pxz, K.T);
    for (size_t i = 0; i < x_dim; i++)
        for (size_t j = 0; j <= i; j++)
            for (size_t k = 0; k < z_dim[z_idx]; k++)
                P(i, j) -= (*Pxz)(i, k) * K(j, k);
    Pxz->release();
}