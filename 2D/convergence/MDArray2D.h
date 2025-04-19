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

    MDArray2D<T> operator-(MDArray2D<T> rhs) const {
      assert((rows_ == rhs.rows_) && (cols_ == rhs.cols_));
      MDArray2D<T> res(rows_, cols_);
      for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
          res(i, j) = (*this)(i,j) - rhs(i, j);
          // printf("%f\n", (*this)(i,j));
          // printf("%f\n", rhs(i,j));
        }
      }
      return res;
    }

    double euclidian_norm() const {
      double res = 0;
      for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
          double val = (*this)(i, j);
          res += val * val;
        }
      }
      // normalize result by num elements
      return res / (rows_ * cols_);
    }
    
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }

    // Get pointer to underlying data.
    T* data_ptr() { return data.data(); }
    const T* data_ptr() const { return data.data(); }
};

#endif // MDARRAY2D_H