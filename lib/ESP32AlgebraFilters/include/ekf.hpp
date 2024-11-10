#ifndef EKF_HPP
#define EKF_HPP

/*
The Ekf class implements the Extended Kalman Filter (EKF), a recursive filter used for estimating the state of a dynamic system from noisy measurements. 
This updated version introduces a more streamlined and flexible approach to EKF implementation compared to the previous version.

Key changes:
- Improved modularity: The class now uses templates and function pointers for greater flexibility in defining state transition and measurement functions.
- Reduced memory usage: Temporary matrices and vectors are allocated dynamically as needed, optimizing memory consumption.
- Simplified interface: The class provides template-based methods for prediction, measurement, and update steps, allowing for easy customization and extension.

functional methods:
- predict(): Propagates the state using the provided state transition function and command input, along with the specified process noise covariance matrix.
- update(): Updates the state based on the measurement, using the predicted measurement, its associated innovation covariance matrix, and optionally, the measurement Jacobian

Next steps:
- Integration with other filtering techniques such as the Unscented Kalman Filter (UKF) for comparison and performance evaluation.
*/

#include "matrix.hpp"
#include "symMatrix.hpp"
using namespace operators;

// template <typename T>
// using Vector_f3 = internal::tmp<Vector<T>> && (*)(const Vector<T> &, const Vector<T> &, const Vector<T> &);
// template <typename T>
// using Vector_f2 = internal::tmp<Vector<T>> && (*)(const Vector<T> &, const Vector<T> &);
// template <typename T>
// using Matrix_f3 = internal::tmp<Matrix<T>> && (*)(const Vector<T> &, const Vector<T> &, const Vector<T> &);
// template <typename T>
// using Matrix_f2 = internal::tmp<Matrix<T>> && (*)(const Vector<T> &, const Vector<T> &);

template <size_t x_dim, size_t u_dim, size_t c_dim = 1, size_t z_num = 1, typename T = float>
class Ekf
{
    private:
// public:
    
    Vector_f3<T> f = nullptr;          // state transition function
    size_t z_dim[z_num] = {0};         // measurement dimensions
    Vector_f2<T> h[z_num] = {nullptr}; // measurement functions
    Matrix_f3<T> Fx = nullptr;         // state transition Jacobian
    Matrix_f3<T> Fu = nullptr;         // state transition Jacobian
    Matrix_f3<T> Fc = nullptr;         // state transition Jacobian
    Matrix_f2<T> H[z_num] = {nullptr}; // measurement Jacobian
    Vector<T> y[z_num];                     // measurement residuals
    float ds[z_num]={1};                     // measurement squared Mahalanobis distance with the estimated state
    float alpha = 1-5.f/(5+1);                       // filter forgetting factor

    Matrix<T> H_val[z_num];                 // measurement Jacobian values
    ldl_matrix<T> S[z_num];                 // innovation covariance
    Vector<T> h_val[z_num];                     // predicted measurement
    Vector<T> prev_X = Vector<T>(x_dim);        // previous state
    Matrix<T> Fx_val = Matrix<T>(x_dim, x_dim); // state transition matrices
    Matrix<T> Fu_val = Matrix<T>(u_dim, x_dim); // state transition matrices
    Matrix<T> Fc_val = Matrix<T>(c_dim, x_dim); // state transition matrices
    rowMajorMatrix<T> K;

    Matrix<T> H_P[z_num];
    inline void finite_diff_Fx(const size_t i, const T eps = 1e-4);
    inline void finite_diff_Fu(const size_t i, const T eps = 1e-4);
    void finite_diff_Fx(){for (size_t i = 0; i < x_dim; i++) finite_diff_Fx(i);};
    void finite_diff_Fu(){for (size_t i = 0; i < u_dim; i++) finite_diff_Fu(i);};

    void finite_diff_H(const size_t z_idx, const size_t i, const T eps = 1e-4);
    void finite_diff_H(const size_t z_idx){for (size_t i = 0; i < x_dim; i++) finite_diff_H(z_idx, i);};
    void compute_S(const size_t z_idx);

    bool initted = false;
public:
    Vector<T> X = Vector<T>(x_dim);              // state
    symMatrix<T> P = symMatrix<T>(x_dim, x_dim); // state covariance
    Vector<T> U = Vector<T>(u_dim);              // command input
    symMatrix<T> Q = symMatrix<T>(u_dim);        // process noise covariance
    Vector<T> C = Vector<T>(c_dim);                  // system parameters
    Ekf();
    Ekf(Vector_f3<T> f) : Ekf() { setPredictionFunction(f); }
    ~Ekf(){};

    void setPredictionFunction(Vector_f3<T> f) { this->f = f; }
    void setJacobianFunction_Fx(Matrix_f3<T> Fx) { this->Fx = Fx; }
    void setJacobianFunction_Fu(Matrix_f3<T> Fu) { this->Fu = Fu; }

    void setMeasurementFunction(Vector_f2<T> h, size_t z_dim, size_t z_idx = 0);
    void setJacobianFunction_H(Matrix_f2<T> H, size_t z_idx = 0) { this->H[z_idx] = H; }

    inline void predict();

    inline void update(const Vector<T> &Z, const symMatrix<T> &R, const size_t z_idx = 0);

    inline float getMahalanobisDistance(const size_t z_idx = 0) const { return ds[z_idx]; }
};

#ifndef EKF_CPP
#include "ekf.cpp"
#endif

#endif // EKF_HPP