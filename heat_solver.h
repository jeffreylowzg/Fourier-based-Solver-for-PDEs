// heat_solver.h
#ifndef HEAT_SOLVER_H
#define HEAT_SOLVER_H

#include <vector>
#include <cassert>
#include <string>

//-------------------------------------------------------------------
// MDArray class for 1D arrays (templated)
//-------------------------------------------------------------------
template<typename T>
class MDArray {
private:
    size_t size_;
    std::vector<T> data;
public:
    MDArray(size_t n) : size_(n), data(n, T()) {}
    
    T& operator[](size_t i) { 
        assert(i < size_);
        return data[i]; 
    }
    const T& operator[](size_t i) const { 
        assert(i < size_);
        return data[i]; 
    }
    size_t size() const { return size_; }
    
    // Helper to get pointer to underlying data (useful for FFTW)
    T* data_ptr() { return data.data(); }
    const T* data_ptr() const { return data.data(); }
};

//-------------------------------------------------------------------
// Function declarations for helper routines (FD and spectral)
//-------------------------------------------------------------------
double laplacian(const MDArray<double>& u, size_t i, double dx);
void RK4_step(MDArray<double>& u, double dt, double dx, double alpha);
void saveSolution(const MDArray<double>& u, const std::string& filename, double dx, double time);

// Spectral step (assumes periodic BC)
// For the heat equation: u_t = alpha u_xx, the Fourier mode û_k evolves as:
//   û_k(t+dt) = û_k(t) * exp(-alpha * k^2 * dt)
// L is the domain length.
void spectral_step(MDArray<double>& u, double dt, double alpha, double L);

#endif // HEAT_SOLVER_H
