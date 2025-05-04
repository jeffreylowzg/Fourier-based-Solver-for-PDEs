/**
 * \ingroup convection_2d
 * \file MDArray2D.h
 * \brief Header for 2D MDArray.
 */

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
    /// @brief constructor
    /// @param rows number of rows
    /// @param cols number of columns
    MDArray2D(size_t rows, size_t cols)
        : rows_(rows), cols_(cols), data(rows * cols, T()) {}

    /// @brief Access element (i,j) with row-major ordering.
    /// @param i row
    /// @param j col
    /// @return reference to element at (i,j)
    T& operator()(size_t i, size_t j) {
        assert(i < rows_ && j < cols_);
        return data[i * cols_ + j];
    }
    /// @brief const access to element
    /// @param i row
    /// @param j col
    /// @return const reference to element
    const T& operator()(size_t i, size_t j) const {
        assert(i < rows_ && j < cols_);
        return data[i * cols_ + j];
    }
    
    /// @brief get num rows
    /// @return num rows
    size_t rows() const { return rows_; }
    /// @brief get num cols
    /// @return num cols
    size_t cols() const { return cols_; }

    /// @brief Get pointer to underlying data.
    /// @return data pointer
    T* data_ptr() { return data.data(); }

    /// @brief const access to data
    /// @return data pointer
    const T* data_ptr() const { return data.data(); }
};

#endif // MDARRAY2D_H