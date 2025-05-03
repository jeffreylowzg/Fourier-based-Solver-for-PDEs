/**
 * \ingroup heat_solver_1d
 * \file heat_solver.cpp
 * \brief this file contains implementations for FD and spectral methods for solving the heat equation
 */
 
/** heat_solver.cpp
 * this file contains implementations for FD and spectral methods for solving
 * the heat equation
 */ 
#include "heat_solver.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <fftw3.h>

/**
 * @brief Compute Laplacian using central differences (Dirichlet BC).
 * 
 * @param u Input array.
 * @param i Index at which to compute the Laplacian.
 * @param dx Grid spacing.
 * @return Approximated Laplacian at index i.
 */
double laplacian(const MDArray<double>& u, size_t i, double dx) {
    if(i == 0 || i == u.size()-1)
        return 0.0;
    return (u[i+1] - 2*u[i] + u[i-1]) / (dx * dx);
}

/**
 * @brief Perform a single RK4 integration step using finite differences.
 * 
 * @param u Array to update.
 * @param dt Time step.
 * @param dx Grid spacing.
 * @param alpha Diffusion coefficient.
 */
void RK4_step(MDArray<double>& u, double dt, double dx, double alpha) {
    size_t n = u.size();
    MDArray<double> k1(n), k2(n), k3(n), k4(n), temp(n);
    
    for(size_t i = 0; i < n; i++) 
        k1[i] = alpha * laplacian(u, i, dx);

    for(size_t i = 0; i < n; i++) 
        temp[i] = u[i] + 0.5 * dt * k1[i];
    for(size_t i = 0; i < n; i++) 
        k2[i] = alpha * laplacian(temp, i, dx);

    for(size_t i = 0; i < n; i++) 
        temp[i] = u[i] + 0.5 * dt * k2[i];
    for(size_t i = 0; i < n; i++) 
        k3[i] = alpha * laplacian(temp, i, dx);

    for(size_t i = 0; i < n; i++) 
        temp[i] = u[i] + dt * k3[i];
    for(size_t i = 0; i < n; i++) 
        k4[i] = alpha * laplacian(temp, i, dx);

    for(size_t i = 0; i < n; i++)
        u[i] += dt / 6.0 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}

/**
 * @brief Save the solution array to a file for visualization.
 * 
 * @param u Solution array.
 * @param filename Output filename.
 * @param dx Grid spacing.
 * @param time Current simulation time.
 */
void saveSolution(const MDArray<double>& u, const std::string& filename, double dx, double time) {
    std::ofstream file("data/" + filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file data/" << filename << " for writing." << std::endl;
        return;
    }
    file << "# Time: " << time << std::endl;
    for(size_t i = 0; i < u.size(); i++){
        double x = i * dx;
        file << x << " " << u[i] << "\n";
    }
    file.close();
}

/**
 * @brief Perform a single RK4 step in Fourier space using the spectral method.
 * 
 * @param u Array to update.
 * @param dt Time step.
 * @param alpha Diffusion coefficient.
 * @param L Domain length.
 */
void spectral_RK4_step(MDArray<double>& u, double dt, double alpha, double L) {
    int n = u.size();
    double* in = fftw_alloc_real(n);
    fftw_complex* out = fftw_alloc_complex(n/2 + 1);

    for (int i = 0; i < n; i++) {
        in[i] = u[i];
    }
    fftw_plan forward = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    for (int k = 0; k < n/2 + 1; k++) {
        double k_val = (2.0 * M_PI * k) / L;
        double lambda = -alpha * k_val * k_val;
        double a = out[k][0];
        double b = out[k][1];

        double k1_a = lambda * a;
        double k1_b = lambda * b;
        double a_k2 = a + 0.5 * dt * k1_a;
        double b_k2 = b + 0.5 * dt * k1_b;
        double k2_a = lambda * a_k2;
        double k2_b = lambda * b_k2;
        double a_k3 = a + 0.5 * dt * k2_a;
        double b_k3 = b + 0.5 * dt * k2_b;
        double k3_a = lambda * a_k3;
        double k3_b = lambda * b_k3;
        double a_k4 = a + dt * k3_a;
        double b_k4 = b + dt * k3_b;
        double k4_a = lambda * a_k4;
        double k4_b = lambda * b_k4;

        out[k][0] = a + dt/6.0 * (k1_a + 2*k2_a + 2*k3_a + k4_a);
        out[k][1] = b + dt/6.0 * (k1_b + 2*k2_b + 2*k3_b + k4_b);
    }

    fftw_plan backward = fftw_plan_dft_c2r_1d(n, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    for (int i = 0; i < n; i++) {
        u[i] = in[i] / n;
    }

    fftw_free(in);
    fftw_free(out);
}

/**
 * @brief Perform a single Backward Euler step in Fourier space using the spectral method.
 * 
 * @param u Array to update.
 * @param dt Time step.
 * @param alpha Diffusion coefficient.
 * @param L Domain length.
 */
void spectral_BE_step(MDArray<double>& u, double dt, double alpha, double L) {
    int n = u.size();
    double* in = fftw_alloc_real(n);
    fftw_complex* out = fftw_alloc_complex(n/2 + 1);

    for (int i = 0; i < n; i++) {
        in[i] = u[i];
    }
    fftw_plan forward = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    for (int k = 0; k < n/2 + 1; k++) {
        double k_val = (2.0 * M_PI * k) / L;
        double factor = 1.0 / (1.0 + alpha * k_val * k_val * dt);
        out[k][0] *= factor;
        out[k][1] *= factor;
    }

    fftw_plan backward = fftw_plan_dft_c2r_1d(n, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    for (int i = 0; i < n; i++) {
        u[i] = in[i] / n;
    }

    fftw_free(in);
    fftw_free(out);
}
