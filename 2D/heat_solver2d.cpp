#include "heat_solver2d.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <fftw3.h>

//-------------------------------------------------------------------
// Source term definitions
//-------------------------------------------------------------------
double source_term(double x, double y, double t, double L) {
    double cx = L / 2.0; // center x (e.g., at origin)
    double cy = L / 2.0; // center y (e.g., at origin)
    double r = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
    double radius = 0.05 * L; // source region radius
    double amplitude = 0.2;   // source amplitude
    double T = 5.0;          // oscillation period
    
    if (r <= radius) {
        // This formulation oscillates between 0 and amplitude.
        // (1 - cos(2*M_PI*t/T)) / 2 ranges from 0 to 1.
        return amplitude * (1 - cos(2 * M_PI * t / T)) / 2.0;
    }
    return 0.0;
}

//-------------------------------------------------------------------
// Compute 2D Laplacian using central differences (Dirichlet BC)
//-------------------------------------------------------------------
double laplacian2d(const MDArray2D<double>& u, size_t i, size_t j, double dx) {
    size_t rows = u.rows();
    size_t cols = u.cols();
    if (i == 0 || i == rows - 1 || j == 0 || j == cols - 1) {
        return 0.0;
    }
    double d2x = (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) / (dx * dx);
    double d2y = (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1)) / (dx * dx);
    return d2x + d2y;
}

//-------------------------------------------------------------------
// RK4 (Explicit) integration step for 2D heat equation (finite-difference)
//-------------------------------------------------------------------
void RK4_step_2d(MDArray2D<double>& u, double dt, double dx, double alpha) {
    size_t rows = u.rows();
    size_t cols = u.cols();
    
    MDArray2D<double> k1(rows, cols), k2(rows, cols), k3(rows, cols), k4(rows, cols), temp(rows, cols);
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            k1(i, j) = alpha * laplacian2d(u, i, j, dx);
        }
    }
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            temp(i, j) = u(i, j) + 0.5 * dt * k1(i, j);
        }
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            k2(i, j) = alpha * laplacian2d(temp, i, j, dx);
        }
    }
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            temp(i, j) = u(i, j) + 0.5 * dt * k2(i, j);
        }
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            k3(i, j) = alpha * laplacian2d(temp, i, j, dx);
        }
    }
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            temp(i, j) = u(i, j) + dt * k3(i, j);
        }
    }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            k4(i, j) = alpha * laplacian2d(temp, i, j, dx);
        }
    }
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            u(i, j) += dt / 6.0 * (k1(i, j) + 2 * k2(i, j) + 2 * k3(i, j) + k4(i, j));
        }
    }
}

//-------------------------------------------------------------------
// Save the 2D solution to a file (suitable for visualization)
//-------------------------------------------------------------------
void saveSolution2D(const MDArray2D<double>& u, const std::string& filename, double dx, double time) {
    std::ofstream file("data/" + filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file data/" << filename << " for writing." << std::endl;
        return;
    }
    file << "# Time: " << time << "\n";
    size_t rows = u.rows();
    size_t cols = u.cols();
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double x = i * dx;
            double y = j * dx;
            file << x << " " << y << " " << u(i, j) << "\n";
        }
        file << "\n";
    }
    file.close();
}

//-------------------------------------------------------------------
// Spectral RK4 step using FFTW (Forward integration in Fourier space)
//-------------------------------------------------------------------
void spectral_RK4_step_2d(MDArray2D<double>& u, double dt, double alpha, double L) {
    int n = u.rows();
    int nc = n;  // n x n grid.
    // Allocate FFTW arrays.
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));
    
    // Copy u into in.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i*nc+j] = u(i, j);
        }
    }
    
    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);
    
    // Loop over Fourier modes.
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq = kx * kx + ky * ky;
            double lambda = -alpha * ksq;
            
            int index = i * (nc/2 + 1) + j;
            double a = out[index][0];
            double b = out[index][1];
            
            // RK4 stages in Fourier space: dU/dt = lambda U.
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
            
            out[index][0] = a + dt / 6.0 * (k1_a + 2*k2_a + 2*k3_a + k4_a);
            out[index][1] = b + dt / 6.0 * (k1_b + 2*k2_b + 2*k3_b + k4_b);
        }
    }
    
    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);
    
    // Normalize and copy back into u.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i*nc+j] / (nc * nc);
        }
    }
    fftw_free(in);
    fftw_free(out);
}

//-------------------------------------------------------------------
// Spectral Backward Euler step using FFTW (Backward integration in Fourier space)
//-------------------------------------------------------------------
void spectral_BE_step_2d(MDArray2D<double>& u, double dt, double alpha, double L) {
    int n = u.rows();
    int nc = n;
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));
    
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i*nc+j] = u(i, j);
        }
    }
    
    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);
    
    // Update each Fourier mode with Backward Euler:
    // (1 - dt * lambda) * U^(n+1) = U^(n)  with lambda = -alpha*k^2, so 1 - dt*lambda = 1 + dt*alpha*k^2.
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq = kx*kx + ky*ky;
            double factor = 1.0 / (1.0 + dt * alpha * ksq);
            int index = i * (nc/2 + 1) + j;
            out[index][0] *= factor;
            out[index][1] *= factor;
        }
    }
    
    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);
    
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i*nc+j] / (nc * nc);
        }
    }
    fftw_free(in);
    fftw_free(out);
}

//-------------------------------------------------------------------
// RK4 step for the 2D heat equation with a source term (finite-difference)
// This computes u_t = alpha * Laplacian(u) + S(x,y,t) with time-dependent S.
//-------------------------------------------------------------------
void RK4_step_2d_source(MDArray2D<double>& u, double dt, double dx, double alpha, double t, double L) {
    size_t rows = u.rows();
    size_t cols = u.cols();
    
    MDArray2D<double> k1(rows, cols), k2(rows, cols), k3(rows, cols), k4(rows, cols), temp(rows, cols);
    
    // Stage 1: f(u, t) = alpha*Laplacian(u) + S(x,y,t)
    for (size_t i = 0; i < rows; i++) {
        double x = i * dx;
        for (size_t j = 0; j < cols; j++) {
            double y = j * dx;
            k1(i, j) = alpha * laplacian2d(u, i, j, dx) + source_term(x, y, t, L);
        }
    }
    
    // Stage 2: f(u + dt/2*k1, t + dt/2)
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            temp(i, j) = u(i, j) + 0.5 * dt * k1(i, j);
        }
    }
    for (size_t i = 0; i < rows; i++) {
        double x = i * dx;
        for (size_t j = 0; j < cols; j++) {
            double y = j * dx;
            k2(i, j) = alpha * laplacian2d(temp, i, j, dx) + source_term(x, y, t + 0.5 * dt, L);
        }
    }
    
    // Stage 3: f(u + dt/2*k2, t + dt/2)
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            temp(i, j) = u(i, j) + 0.5 * dt * k2(i, j);
        }
    }
    for (size_t i = 0; i < rows; i++) {
        double x = i * dx;
        for (size_t j = 0; j < cols; j++) {
            double y = j * dx;
            k3(i, j) = alpha * laplacian2d(temp, i, j, dx) + source_term(x, y, t + 0.5 * dt, L);
        }
    }
    
    // Stage 4: f(u + dt*k3, t + dt)
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            temp(i, j) = u(i, j) + dt * k3(i, j);
        }
    }
    for (size_t i = 0; i < rows; i++) {
        double x = i * dx;
        for (size_t j = 0; j < cols; j++) {
            double y = j * dx;
            k4(i, j) = alpha * laplacian2d(temp, i, j, dx) + source_term(x, y, t + dt, L);
        }
    }
    
    // Update the solution.
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            u(i, j) += dt / 6.0 * (k1(i, j) + 2*k2(i, j) + 2*k3(i, j) + k4(i, j));
        }
    }
}

//-------------------------------------------------------------------
// Spectral RK4 step with source using FFTW.
// This function integrates u_t = alpha * Laplacian(u) + S(x,y,t)
// by first applying the spectral RK4 update for the diffusion part
// and then adding the source contribution computed via RK4 in physical space.
//-------------------------------------------------------------------
void spectral_RK4_step_2d_source(MDArray2D<double>& u, double dt, double alpha, double L, double t) {
    int n = u.rows();
    int nc = n;  // n x n grid.
    // Allocate FFTW arrays.
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));
    
    // Copy u into in.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i*nc+j] = u(i, j);
        }
    }
    
    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);
    
    // Perform RK4 in Fourier space for the diffusion term.
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq = kx * kx + ky * ky;
            double lambda = -alpha * ksq;
            
            int index = i * (nc/2 + 1) + j;
            double a = out[index][0];
            double b = out[index][1];
            
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
            
            out[index][0] = a + dt / 6.0 * (k1_a + 2*k2_a + 2*k3_a + k4_a);
            out[index][1] = b + dt / 6.0 * (k1_b + 2*k2_b + 2*k3_b + k4_b);
        }
    }
    
    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);
    
    // Normalize and copy into a temporary physical-space array.
    MDArray2D<double> u_diff(nc, nc);
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u_diff(i, j) = in[i*nc+j] / (nc * nc);
        }
    }
    fftw_free(in);
    fftw_free(out);
    
    // Compute the source contribution in physical space using RK4 integration.
    double dx = L / (n - 1);
    for (size_t i = 0; i < u_diff.rows(); i++) {
        double x = i * dx;
        for (size_t j = 0; j < u_diff.cols(); j++) {
            double y = j * dx;
            double k1 = source_term(x, y, t, L);
            double k2 = source_term(x, y, t + 0.5 * dt, L);
            double k3 = source_term(x, y, t + 0.5 * dt, L);
            double k4 = source_term(x, y, t + dt, L);
            double source_update = dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
            u_diff(i, j) += source_update;
        }
    }
    
    // Copy the updated result back into u.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = u_diff(i, j);
        }
    }
}

//-------------------------------------------------------------------
// Spectral Backward Euler step with source using FFTW.
// This function integrates u_t = alpha * Laplacian(u) + S(x,y,t)
// by first applying the spectral Backward Euler update for the diffusion part
// and then adding the source contribution computed via RK4 in physical space.
//-------------------------------------------------------------------
void spectral_BE_step_2d_source(MDArray2D<double>& u, double dt, double alpha, double L, double t) {
    int n = u.rows();
    int nc = n;
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));
    
    // Copy u into in.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i*nc+j] = u(i, j);
        }
    }
    
    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);
    
    // Apply Backward Euler update for diffusion in Fourier space.
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq = kx * kx + ky * ky;
            double factor = 1.0 / (1.0 + dt * alpha * ksq);
            int index = i * (nc/2 + 1) + j;
            out[index][0] *= factor;
            out[index][1] *= factor;
        }
    }
    
    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);
    
    MDArray2D<double> u_diff(nc, nc);
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u_diff(i, j) = in[i*nc+j] / (nc * nc);
        }
    }
    fftw_free(in);
    fftw_free(out);
    
    // Compute the source contribution in physical space using RK4 integration.
    double dx = L / (n - 1);
    for (size_t i = 0; i < u_diff.rows(); i++) {
        double x = i * dx;
        for (size_t j = 0; j < u_diff.cols(); j++) {
            double y = j * dx;
            double k1 = source_term(x, y, t, L);
            double k2 = source_term(x, y, t + 0.5 * dt, L);
            double k3 = source_term(x, y, t + 0.5 * dt, L);
            double k4 = source_term(x, y, t + dt, L);
            double source_update = dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
            u_diff(i, j) += source_update;
        }
    }
    
    // Copy the updated result back into u.
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = u_diff(i, j);
        }
    }
}
