#ifndef MDARRAY2D_H
#define MDARRAY2D_H

#include <vector>
#include <cassert>
#include <cstddef>

// Templated multi-dimensional array class for 2D arrays.
template<typename T>
class MDArray2D {
private:
    size_t rows_, cols_;
    std::vector<T> data;
public:
    MDArray2D(size_t rows, size_t cols)
        : rows_(rows), cols_(cols), data(rows * cols, T()) {}

    // Access element (i,j) with row-major ordering.
    T& operator()(size_t i, size_t j) {
        assert(i < rows_ && j < cols_);
        return data[i * cols_ + j];
    }
    const T& operator()(size_t i, size_t j) const {
        assert(i < rows_ && j < cols_);
        return data[i * cols_ + j];
    }
    
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }

    // Get pointer to underlying data.
    T* data_ptr() { return data.data(); }
    const T* data_ptr() const { return data.data(); }
};

#endif // MDARRAY2D_H