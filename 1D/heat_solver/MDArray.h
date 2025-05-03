/**
 * \ingroup heat_solver_1d
 * \file MDArray.h
 * \brief Header file to handle multi-dimensional array.
 */

/** MDArray.h
 * 
 * header file for 1D MDArray
 */

#ifndef MDARRAY_H
#define MDARRAY_H

#include <vector>
#include <cassert>
#include <cstddef> // for size_t

/**
 * @brief Templated 1D array class with basic access and utility functions.
 * 
 * @tparam T Type of the array elements.
 */
template<typename T>
class MDArray {
private:
    size_t size_;
    std::vector<T> data;
public:
    MDArray(size_t n) : size_(n), data(n, T()) {}

    /**
     * @brief Access element at index i (non-const).
     * 
     * @param i Index.
     * @return Reference to element at index i.
     */
    T& operator[](size_t i) {
        assert(i < size_);
        return data[i];
    }

    /**
     * @brief Access element at index i (const).
     * 
     * @param i Index.
     * @return Const reference to element at index i.
     */
    const T& operator[](size_t i) const {
        assert(i < size_);
        return data[i];
    }

    /**
     * @brief Get the size of the array.
     * 
     * @return Number of elements.
     */
    size_t size() const { return size_; }

    /**
     * @brief Get pointer to underlying data array (non-const).
     * 
     * @return Pointer to data.
     */
    T* data_ptr() { return data.data(); }

    /**
     * @brief Get pointer to underlying data array (const).
     * 
     * @return Const pointer to data.
     */
    const T* data_ptr() const { return data.data(); }
};

#endif // MDARRAY_H
