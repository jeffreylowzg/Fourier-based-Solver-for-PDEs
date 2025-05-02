/*
 * heat_solver2d.cpp
 * Adapted to support selectable source-term profiles.
 * Original structure and style preserved.
 *
 * This code uses the pseudo-spectral method to solve with source terms
 * The source terms you can use are:
 *  - gaussian pulse
 *  - standing wave
 *  - traveling gaussian
 *  - pulsed ring
 * 
 * You can use spectral backward Euler or RK4
 * Or you can use finite difference RK4
 * 
 */

#include "heat_solver2d.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <string>
#include <fftw3.h>

//-------------------------------------------------------------------
// Source term definitions
//-------------------------------------------------------------------
double source_term(double x, double y, double t, double L, const std::string& source_type) {
    double cx = 0.5 * L;
    double cy = 0.5 * L;

    if (source_type == "gaussian_pulse") {
        double A     = 1.0;
        double sigma = 0.1 * L; //controls the diameter
        double tau   = 1.0;
        double dx2   = (x - cx)*(x - cx);
        double dy2   = (y - cy)*(y - cy);
        double spatial  = std::exp(-(dx2 + dy2)/(2.0*sigma*sigma));
        double temporal = 1.0 - std::exp(-t/tau);
        return A * spatial * temporal;
    }
    else if (source_type == "standing_wave") {
        double A      = 1.0;
        double period = 5.0;
        double omega  = 2.0 * M_PI / period;
        int mx = 4; 
        int my = 4;   // control size of wave
        return A
        * sin(mx*M_PI*x/L)
        * sin(my*M_PI*y/L)
        * cos(omega*t);

    }
    else if (source_type == "traveling_gaussian") {
        double A     = 1.0;
        double v     = 0.2 * L;
        double sigma = 0.05 * L;
        double x0    = 0.1 * L;
        double xc    = x0 + v * t;
        double dx2   = (x - xc)*(x - xc);
        double dy2   = (y - cy)*(y - cy);
        return A * std::exp(-(dx2 + dy2)/(2.0*sigma*sigma));
    }
    else if (source_type == "pulsed_ring") {
        double r      = std::sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy));
        double r1     = 0.2 * L;
        double r2     = 0.3 * L;
        double t0     = 1.0;
        double dt_pulse = 0.5;
        if (r >= r1 && r <= r2 && t >= t0 && t <= t0 + dt_pulse)
            return 1.0;
        else
            return 0.0;
    }
    // none or unrecognized type
    return 0.0;
}

//-------------------------------------------------------------------
// Compute 2D Laplacian using central differences (Dirichlet BC)
//-------------------------------------------------------------------
double laplacian2d(const MDArray2D<double>& u, size_t i, size_t j, double dx) {
    size_t rows = u.rows();
    size_t cols = u.cols();
    if (i == 0 || i == rows-1 || j == 0 || j == cols-1) {
        return 0.0;
    }
    double d2x = (u(i+1,j) - 2.0*u(i,j) + u(i-1,j)) / (dx*dx);
    double d2y = (u(i,j+1) - 2.0*u(i,j) + u(i,j-1)) / (dx*dx);
    return d2x + d2y;
}

//-------------------------------------------------------------------
// RK4 (Explicit) integration step for 2D heat equation (finite-difference)
//-------------------------------------------------------------------
void RK4_step_2d(MDArray2D<double>& u, double dt, double dx, double alpha) {
    size_t rows = u.rows();
    size_t cols = u.cols();
    MDArray2D<double> k1(rows,cols), k2(rows,cols),
                     k3(rows,cols), k4(rows,cols),
                     temp(rows,cols);

    // Stage 1
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            k1(i,j) = alpha * laplacian2d(u, i, j, dx);

    // Stage 2
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            temp(i,j) = u(i,j) + 0.5 * dt * k1(i,j);

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            k2(i,j) = alpha * laplacian2d(temp, i, j, dx);

    // Stage 3
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            temp(i,j) = u(i,j) + 0.5 * dt * k2(i,j);

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            k3(i,j) = alpha * laplacian2d(temp, i, j, dx);

    // Stage 4
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            temp(i,j) = u(i,j) + dt * k3(i,j);

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            k4(i,j) = alpha * laplacian2d(temp, i, j, dx);

    // Combine
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            u(i,j) += dt / 6.0 * (k1(i,j) + 2.0*k2(i,j) + 2.0*k3(i,j) + k4(i,j));
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
            file << x << " " << y << " " << u(i,j) << "\n";
        }
        file << "\n";
    }
    file.close();
}

//-------------------------------------------------------------------
// Spectral RK4 step using FFTW (Forward integration in Fourier space)
//-------------------------------------------------------------------
void spectral_RK4_step_2d(MDArray2D<double>& u, double dt, double alpha, double L) {
    int n  = u.rows();
    int nc = n;
    double* in  = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));

    // Copy u into in
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < nc; j++)
            in[i*nc + j] = u(i,j);

    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    // RK4 in Fourier space
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq   = kx*kx + ky*ky;
            double lambda = -alpha * ksq;
            int index     = i * (nc/2 + 1) + j;
            double a = out[index][0];
            double b = out[index][1];

            double k1_a = lambda * a;
            double k1_b = lambda * b;
            double a2   = a + 0.5 * dt * k1_a;
            double b2   = b + 0.5 * dt * k1_b;
            double k2_a = lambda * a2;
            double k2_b = lambda * b2;
            double a3   = a + 0.5 * dt * k2_a;
            double b3   = b + 0.5 * dt * k2_b;
            double k3_a = lambda * a3;
            double k3_b = lambda * b3;
            double a4   = a + dt * k3_a;
            double b4   = b + dt * k3_b;
            double k4_a = lambda * a4;
            double k4_b = lambda * b4;

            out[index][0] = a + dt / 6.0 * (k1_a + 2*k2_a + 2*k3_a + k4_a);
            out[index][1] = b + dt / 6.0 * (k1_b + 2*k2_b + 2*k3_b + k4_b);
        }
    }

    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    // Normalize and copy back into u
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < nc; j++)
            u(i,j) = in[i*nc + j] / (nc * nc);

    fftw_free(in);
    fftw_free(out);
}

//-------------------------------------------------------------------
// Spectral Backward Euler step using FFTW (Backward integration in Fourier space)
//-------------------------------------------------------------------
void spectral_BE_step_2d(MDArray2D<double>& u, double dt, double alpha, double L) {
    int n  = u.rows();
    int nc = n;
    double* in  = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));

    // Copy u into in
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < nc; j++)
            in[i*nc + j] = u(i,j);

    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    // Backward Euler in Fourier space
    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx     = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky     = 2.0 * M_PI * k_j / L;
            double ksq    = kx*kx + ky*ky;
            double factor = 1.0 / (1.0 + dt * alpha * ksq);
            int index     = i * (nc/2 + 1) + j;
            out[index][0] *= factor;
            out[index][1] *= factor;
        }
    }

    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    // Normalize and copy back into u
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < nc; j++)
            u(i,j) = in[i*nc + j] / (nc * nc);

    fftw_free(in);
    fftw_free(out);
}

//-------------------------------------------------------------------
// RK4 step for the 2D heat equation with a source term (finite-difference)
//    u_t = alpha * Laplacian(u) + S(x,y,t)
//-------------------------------------------------------------------
void RK4_step_2d_source(MDArray2D<double>& u, double dt, double dx,
                        double alpha, double t, double L,
                        const std::string& source_type) {
    size_t rows = u.rows();
    size_t cols = u.cols();
    MDArray2D<double> k1(rows,cols), k2(rows,cols),
                     k3(rows,cols), k4(rows,cols),
                     temp(rows,cols);

    // Stage 1
    for (size_t i = 0; i < rows; i++) {
        double x = i * dx;
        for (size_t j = 0; j < cols; j++) {
            double y = j * dx;
            k1(i,j) = alpha * laplacian2d(u,i,j,dx)
                    + source_term(x,y,t,   L, source_type);
        }
    }
    // Stage 2
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            temp(i,j) = u(i,j) + 0.5 * dt * k1(i,j);

    for (size_t i = 0; i < rows; i++) {
        double x = i * dx;
        for (size_t j = 0; j < cols; j++) {
            double y = j * dx;
            k2(i,j) = alpha * laplacian2d(temp,i,j,dx)
                    + source_term(x,y,t+0.5*dt,L, source_type);
        }
    }
    // Stage 3
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            temp(i,j) = u(i,j) + 0.5 * dt * k2(i,j);

    for (size_t i = 0; i < rows; i++) {
        double x = i * dx;
        for (size_t j = 0; j < cols; j++) {
            double y = j * dx;
            k3(i,j) = alpha * laplacian2d(temp,i,j,dx)
                    + source_term(x,y,t+0.5*dt,L, source_type);
        }
    }
    // Stage 4
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            temp(i,j) = u(i,j) + dt * k3(i,j);

    for (size_t i = 0; i < rows; i++) {
        double x = i * dx;
        for (size_t j = 0; j < cols; j++) {
            double y = j * dx;
            k4(i,j) = alpha * laplacian2d(temp,i,j,dx)
                    + source_term(x,y,t+dt,  L, source_type);
        }
    }
    // Combine
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            u(i,j) += dt / 6.0 *
                     (k1(i,j) + 2.0*k2(i,j) + 2.0*k3(i,j) + k4(i,j));
}

//-------------------------------------------------------------------
// Spectral RK4 step with source using FFTW.
//-------------------------------------------------------------------
void spectral_RK4_step_2d_source(MDArray2D<double>& u,
                                 double dt, double alpha,
                                 double L, double t,
                                 const std::string& source_type) {
    // diffusion-only
    spectral_RK4_step_2d(u, dt, alpha, L);

    // source in physical space
    int n = u.rows();
    double dx = L / (n - 1);
    MDArray2D<double> u_diff = u;
    for (int i = 0; i < n; i++) {
        double x = i * dx;
        for (int j = 0; j < n; j++) {
            double y = j * dx;
            double s1 = source_term(x,y,      t,       L, source_type);
            double s2 = source_term(x,y,      t+0.5*dt,L, source_type);
            double s3 = source_term(x,y,      t+0.5*dt,L, source_type);
            double s4 = source_term(x,y,      t+dt,    L, source_type);
            double ds = dt / 6.0 * (s1 + 2.0*s2 + 2.0*s3 + s4);
            u_diff(i,j) += ds;
        }
    }
    u = u_diff;
}

//-------------------------------------------------------------------
// Spectral Backward Euler step with source using FFTW.
//-------------------------------------------------------------------
void spectral_BE_step_2d_source(MDArray2D<double>& u,
                                double dt, double alpha,
                                double L, double t,
                                const std::string& source_type) {
    // diffusion-only
    spectral_BE_step_2d(u, dt, alpha, L);

    // source in physical space
    int n = u.rows();
    double dx = L / (n - 1);
    MDArray2D<double> u_diff = u;
    for (int i = 0; i < n; i++) {
        double x = i * dx;
        for (int j = 0; j < n; j++) {
            double y = j * dx;
            double s1 = source_term(x,y,      t,       L, source_type);
            double s2 = source_term(x,y,      t+0.5*dt,L, source_type);
            double s3 = source_term(x,y,      t+0.5*dt,L, source_type);
            double s4 = source_term(x,y,      t+dt,    L, source_type);
            double ds = dt / 6.0 * (s1 + 2.0*s2 + 2.0*s3 + s4);
            u_diff(i,j) += ds;
        }
    }
    u = u_diff;
}
