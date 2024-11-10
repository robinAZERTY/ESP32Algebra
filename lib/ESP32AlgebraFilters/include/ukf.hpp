#ifndef UKF_HPP
#define UKF_HPP
#include "matrix.hpp"
#include "symMatrix.hpp"
using namespace operators;
#ifndef ARDUINO
#include <math.h>
#endif

template <size_t x_dim, size_t u_dim, size_t c_dim = 1, size_t z_num = 1, typename T = float>
class Ukf
{
    // private:
public:
    Matrix<T> K;
    
    Vector_f3<T> f = nullptr;          // state transition function
    Matrix<T> X_sp;                     // state sigma points
    size_t z_dim[z_num] = {0};         // measurement dimensions
    Vector_f2<T> h[z_num] = {nullptr}; // measurement functions

    ldl_matrix<T> S[z_num];                 // innovation covariance
    Matrix<T> h_sp[z_num];                     // predicted measurement sigma points
    Vector<T> h_val[z_num];                    // predicted measurement
    void compute_S(const size_t z_idx);
    const T sqrt_x_dim_u_dim = sqrt(x_dim + u_dim);
    void compute_state_sigma_points();
    void compute_measurement_sigma_points(const size_t z_idx);
public:
    Vector<T> X = Vector<T>(x_dim);              // state
    ldl_matrix<T> P = ldl_matrix<T>(x_dim); // state covariance
    Vector<T> U = Vector<T>(u_dim);              // command input
    ldl_matrix<T> Q = ldl_matrix<T>(u_dim);        // process noise covariance
    Vector<T> C = Vector<T>(c_dim);                  // system parameters

    Ukf();
    Ukf(Vector_f3<T> f) : Ukf() { setPredictionFunction(f); }
    ~Ukf(){};

    void setPredictionFunction(Vector_f3<T> f) { this->f = f; }
    void setMeasurementFunction(Vector_f2<T> h, size_t z_dim, size_t z_idx = 0);
    void predict();
    void update(const Vector<T> &Z, const symMatrix<T> &R, const size_t z_idx = 0);
};

#ifndef UKF_CPP
#include "ukf.cpp"
#endif

#endif // UKF_HPP