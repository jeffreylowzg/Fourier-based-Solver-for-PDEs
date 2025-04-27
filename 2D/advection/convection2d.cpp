#include "convection2d.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <fftw3.h>

//-------------------------------------------------------------------
// Spectral RK4 step for the 2D convection–diffusion equation
//-------------------------------------------------------------------
void spectral_RK4_step_2d_convection(MDArray2D<double>& u,
                                     double dt, double D,
                                     double vx, double vy,
                                     double L)
{
    int n = static_cast<int>(u.rows());
    int nc = n;
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i * nc + j] = u(i, j);
        }
    }

    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq = kx * kx + ky * ky;
            double lambda_r = -D * ksq;
            double lambda_i = -(kx * vx + ky * vy);

            int index = i * (nc/2 + 1) + j;
            double a = out[index][0];
            double b = out[index][1];

            if (ksq == 0.0 && lambda_i == 0.0)
                continue;

            double k1_a = lambda_r * a - lambda_i * b;
            double k1_b = lambda_i * a + lambda_r * b;

            double a_k2 = a + 0.5 * dt * k1_a;
            double b_k2 = b + 0.5 * dt * k1_b;
            double k2_a = lambda_r * a_k2 - lambda_i * b_k2;
            double k2_b = lambda_i * a_k2 + lambda_r * b_k2;

            double a_k3 = a + 0.5 * dt * k2_a;
            double b_k3 = b + 0.5 * dt * k2_b;
            double k3_a = lambda_r * a_k3 - lambda_i * b_k3;
            double k3_b = lambda_i * a_k3 + lambda_r * b_k3;

            double a_k4 = a + dt * k3_a;
            double b_k4 = b + dt * k3_b;
            double k4_a = lambda_r * a_k4 - lambda_i * b_k4;
            double k4_b = lambda_i * a_k4 + lambda_r * b_k4;

            out[index][0] = a + dt / 6.0 * (k1_a + 2*k2_a + 2*k3_a + k4_a);
            out[index][1] = b + dt / 6.0 * (k1_b + 2*k2_b + 2*k3_b + k4_b);
        }
    }

    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    double invN2 = 1.0 / (nc * nc);
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i * nc + j] * invN2;
        }
    }

    fftw_free(in);
    fftw_free(out);
}

//-------------------------------------------------------------------
// Spectral Backward Euler step for the 2D convection–diffusion equation
//-------------------------------------------------------------------
void spectral_BE_step_2d_convection(MDArray2D<double>& u,
                                    double dt, double D,
                                    double vx, double vy,
                                    double L)
{
    int n = static_cast<int>(u.rows());
    int nc = n;
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i * nc + j] = u(i, j);
        }
    }

    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    for (int i = 0; i < nc; i++) {
        int k_i = (i <= nc/2) ? i : i - nc;
        double kx = 2.0 * M_PI * k_i / L;
        for (int j = 0; j < nc/2 + 1; j++) {
            int k_j = j;
            double ky = 2.0 * M_PI * k_j / L;
            double ksq = kx * kx + ky * ky;
            double lambda_r = -D * ksq;
            double lambda_i = -(kx * vx + ky * vy);

            int index = i * (nc/2 + 1) + j;
            double a = out[index][0];
            double b = out[index][1];

            if (ksq == 0.0 && lambda_i == 0.0)
                continue;

            double denom_r = 1.0 - dt * lambda_r;
            double denom_i = -dt * lambda_i;
            double denom_mag2 = denom_r * denom_r + denom_i * denom_i;

            double new_a = (a * denom_r + b * denom_i) / denom_mag2;
            double new_b = (b * denom_r - a * denom_i) / denom_mag2;

            out[index][0] = new_a;
            out[index][1] = new_b;
        }
    }

    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    double invN2 = 1.0 / (nc * nc);
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i * nc + j] * invN2;
        }
    }

    fftw_free(in);
    fftw_free(out);
}

//-------------------------------------------------------------------
// FULL Spectral RK4 step for the 2D convection–diffusion equation
//-------------------------------------------------------------------
void full_spectral_RK4_step_2d_convection(MDArray2D<double>& u,
                                          double dt, double D,
                                          double vx, double vy,
                                          double L)
{
    int n = static_cast<int>(u.rows());
    int nc = n;
    double* in = fftw_alloc_real(nc * nc);
    fftw_complex* out = fftw_alloc_complex(nc * (nc/2 + 1));

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            in[i * nc + j] = u(i, j);
        }
    }

    fftw_plan forward = fftw_plan_dft_r2c_2d(nc, nc, in, out, FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    std::vector<std::pair<double, double>> k1(nc*(nc/2+1));
    std::vector<std::pair<double, double>> k2(nc*(nc/2+1));
    std::vector<std::pair<double, double>> k3(nc*(nc/2+1));
    std::vector<std::pair<double, double>> k4(nc*(nc/2+1));

    auto compute_lambda = [&](int i, int j, double& lambda_r, double& lambda_i){
        int k_i = (i <= nc/2) ? i : i - nc;
        int k_j = j;
        double kx = 2.0 * M_PI * k_i / L;
        double ky = 2.0 * M_PI * k_j / L;
        double ksq = kx*kx + ky*ky;
        lambda_r = -D * ksq;
        lambda_i = -(kx * vx + ky * vy);
    };

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= nc/2; j++) {
            double lambda_r, lambda_i;
            compute_lambda(i, j, lambda_r, lambda_i);

            int index = i*(nc/2+1)+j;
            double a = out[index][0];
            double b = out[index][1];

            k1[index].first = lambda_r * a - lambda_i * b;
            k1[index].second = lambda_i * a + lambda_r * b;
        }
    }

    fftw_complex* temp = fftw_alloc_complex(nc*(nc/2+1));

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= nc/2; j++) {
            int index = i*(nc/2+1)+j;
            temp[index][0] = out[index][0] + 0.5 * dt * k1[index].first;
            temp[index][1] = out[index][1] + 0.5 * dt * k1[index].second;
        }
    }

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= nc/2; j++) {
            double lambda_r, lambda_i;
            compute_lambda(i, j, lambda_r, lambda_i);

            int index = i*(nc/2+1)+j;
            double a = temp[index][0];
            double b = temp[index][1];

            k2[index].first = lambda_r * a - lambda_i * b;
            k2[index].second = lambda_i * a + lambda_r * b;
        }
    }

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= nc/2; j++) {
            int index = i*(nc/2+1)+j;
            temp[index][0] = out[index][0] + 0.5 * dt * k2[index].first;
            temp[index][1] = out[index][1] + 0.5 * dt * k2[index].second;
        }
    }

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= nc/2; j++) {
            double lambda_r, lambda_i;
            compute_lambda(i, j, lambda_r, lambda_i);

            int index = i*(nc/2+1)+j;
            double a = temp[index][0];
            double b = temp[index][1];

            k3[index].first = lambda_r * a - lambda_i * b;
            k3[index].second = lambda_i * a + lambda_r * b;
        }
    }

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= nc/2; j++) {
            int index = i*(nc/2+1)+j;
            temp[index][0] = out[index][0] + dt * k3[index].first;
            temp[index][1] = out[index][1] + dt * k3[index].second;
        }
    }

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= nc/2; j++) {
            double lambda_r, lambda_i;
            compute_lambda(i, j, lambda_r, lambda_i);

            int index = i*(nc/2+1)+j;
            double a = temp[index][0];
            double b = temp[index][1];

            k4[index].first = lambda_r * a - lambda_i * b;
            k4[index].second = lambda_i * a + lambda_r * b;
        }
    }

    for (int i = 0; i < nc; i++) {
        for (int j = 0; j <= nc/2; j++) {
            int index = i*(nc/2+1)+j;
            out[index][0] += dt/6.0 * (k1[index].first + 2*k2[index].first + 2*k3[index].first + k4[index].first);
            out[index][1] += dt/6.0 * (k1[index].second + 2*k2[index].second + 2*k3[index].second + k4[index].second);
        }
    }

    fftw_plan backward = fftw_plan_dft_c2r_2d(nc, nc, out, in, FFTW_ESTIMATE);
    fftw_execute(backward);
    fftw_destroy_plan(backward);

    double invN2 = 1.0 / (nc * nc);
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < nc; j++) {
            u(i, j) = in[i * nc + j] * invN2;
        }
    }

    fftw_free(in);
    fftw_free(out);
    fftw_free(temp);
}

//-------------------------------------------------------------------
// Save the 2D solution to a file
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
