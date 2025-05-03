/**
 * \ingroup convergence_1d
 * \file MDArray.h
 * \brief Header file to handle multi-dimensional array.
 */

#ifndef MDARRAY_H
#define MDARRAY_H

#include <vector>
#include <cassert>
#include <cstddef> // for size_t

/**
 * @brief Templated 1D multi-dimensional array class.
 * 
 * @tparam T Type of the array elements.
 */
template<typename T>
class MDArray {
private:
    size_t size_;
    std::vector<T> data;
public:
    /**
     * @brief Constructor to initialize array with size n.
     * 
     * @param n Number of elements.
     */
    MDArray(size_t n) : size_(n), data(n, T()) {}

    /**
     * @brief Access element at index i (non-const).
     * 
     * @param i Index.
     * @return Reference to element.
     */
    T& operator[](size_t i) {
        assert(i < size_);
        return data[i];
    }

    /**
     * @brief Access element at index i (const).
     * 
     * @param i Index.
     * @return Const reference to element.
     */
    const T& operator[](size_t i) const {
        assert(i < size_);
        return data[i];
    }

    /**
     * @brief Get the number of elements in the array.
     * 
     * @return Number of elements.
     */
    size_t size() const { return size_; }

    /**
     * @brief Subtract another MDArray and return the result.
     * 
     * @param rhs Right-hand side array.
     * @return Resulting MDArray after subtraction.
     */
    MDArray<T> operator-(MDArray<T> rhs) const {
        assert(size_ == rhs.size_);
        MDArray<T> res(size_);
        for (size_t i = 0; i < size_; i++) {
            res[i] = (*this)[i] - rhs[i];
        }
        return res;
    }

    /**
     * @brief Calculate the Euclidean norm (L2 norm) of the array.
     * 
     * @return Normalized Euclidean norm.
     */
    double euclidian_norm() const {
        double res = 0;
        for (size_t i = 0; i < size_; i++) {
            double val = (*this)[i];
            res += val * val;
        }
        return res / (size_);
    }

    /**
     * @brief Get a pointer to the underlying data (non-const).
     * 
     * @return Pointer to data.
     */
    T* data_ptr() { return data.data(); }

    /**
     * @brief Get a pointer to the underlying data (const).
     * 
     * @return Const pointer to data.
     */
    const T* data_ptr() const { return data.data(); }
};

#endif // MDARRAY_H
