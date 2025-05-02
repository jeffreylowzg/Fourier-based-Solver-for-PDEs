/** MDArray.h
 * 
 * header file for 1D MDArray
 */

#ifndef MDARRAY_H
#define MDARRAY_H

#include <vector>
#include <cassert>
#include <cstddef> // for size_t

// Templated multi-dimensional array class for 1D arrays
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

#endif // MDARRAY_H
