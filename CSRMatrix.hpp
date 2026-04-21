#ifndef CSR_MATRIX_HPP
#define CSR_MATRIX_HPP

#include <vector>
#include <exception>

namespace sjtu {

class size_mismatch : public std::exception {
public:
    const char *what() const noexcept override {
        return "Size mismatch";
    }
};

class invalid_index : public std::exception {
public:
    const char *what() const noexcept override {
        return "Index out of range";
    }
};

// CSR (Compressed Sparse Row) matrix implementation
// Only std::vector is used as requested.
template <typename T>
class CSRMatrix {

private:
    size_t nrows_{};
    size_t ncols_{};
    std::vector<size_t> indptr_;
    std::vector<size_t> indices_;
    std::vector<T> data_;

    // find position of column j within row i using binary search
    // returns index into indices_/data_ or (size_t)-1 if not found
    size_t find_in_row(size_t i, size_t j) const {
        size_t l = indptr_[i];
        size_t r = indptr_[i + 1];
        // Binary search since we maintain sorted indices per row
        while (l < r) {
            size_t mid = l + (r - l) / 2;
            size_t col = indices_[mid];
            if (col == j) return mid;
            if (col < j) l = mid + 1; else r = mid;
        }
        return static_cast<size_t>(-1);
    }

    // find insertion position in row i to keep indices sorted
    size_t lower_bound_in_row(size_t i, size_t j) const {
        size_t l = indptr_[i];
        size_t r = indptr_[i + 1];
        while (l < r) {
            size_t mid = l + (r - l) / 2;
            if (indices_[mid] < j) l = mid + 1; else r = mid;
        }
        return l; // first position with col >= j
    }

public:
    CSRMatrix &operator=(const CSRMatrix &other) = delete;
    CSRMatrix &operator=(CSRMatrix &&other) = delete;

    // Initialize an empty CSR matrix with n rows and m columns
    CSRMatrix(size_t n, size_t m) : nrows_(n), ncols_(m), indptr_(n + 1, 0) {}

    // Initialize CSR matrix from existing CSR format data, validate sizes
    CSRMatrix(size_t n, size_t m, size_t count,
              const std::vector<size_t> &indptr,
              const std::vector<size_t> &indices,
              const std::vector<T> &data)
        : nrows_(n), ncols_(m), indptr_(indptr), indices_(indices), data_(data) {
        if (indptr_.size() != nrows_ + 1) {
            throw size_mismatch();
        }
        if (indices_.size() != count || data_.size() != count) {
            throw size_mismatch();
        }
        // indptr must be non-decreasing and last equals count
        if (indptr_.empty() || indptr_[0] != 0 || indptr_.back() != count) {
            throw size_mismatch();
        }
        for (size_t i = 1; i < indptr_.size(); ++i) {
            if (indptr_[i] < indptr_[i - 1] || indptr_[i] > count) throw size_mismatch();
        }
        // indices must be within column range
        for (size_t k = 0; k < indices_.size(); ++k) {
            if (indices_[k] >= ncols_) throw invalid_index();
        }
        // ensure column indices are strictly increasing within each row
        for (size_t i = 0; i < nrows_; ++i) {
            size_t s = indptr_[i], e = indptr_[i + 1];
            for (size_t p = s; p + 1 < e; ++p) {
                if (!(indices_[p] < indices_[p + 1])) throw size_mismatch();
            }
        }
    }

    // Copy constructor
    CSRMatrix(const CSRMatrix &other) = default;

    // Move constructor
    CSRMatrix(CSRMatrix &&other) = default;

    // Convert dense matrix representation to CSR format
    CSRMatrix(size_t n, size_t m, const std::vector<std::vector<T>> &dense)
        : nrows_(n), ncols_(m), indptr_(n + 1, 0) {
        if (dense.size() != nrows_) throw size_mismatch();
        indices_.clear();
        data_.clear();
        for (size_t i = 0; i < nrows_; ++i) {
            const auto &row = dense[i];
            if (row.size() != ncols_) throw size_mismatch();
            for (size_t j = 0; j < ncols_; ++j) {
                if (!(row[j] == T{})) {
                    indices_.push_back(j);
                    data_.push_back(row[j]);
                }
            }
            indptr_[i + 1] = indices_.size();
        }
    }

    ~CSRMatrix() = default;

    // Return the number of rows
    size_t getRowSize() const { return nrows_; }

    // Return the number of columns
    size_t getColSize() const { return ncols_; }

    // Return the count of non-zero elements
    size_t getNonZeroCount() const { return data_.size(); }

    // Retrieve element at position (i,j)
    T get(size_t i, size_t j) const {
        if (i >= nrows_ || j >= ncols_) throw invalid_index();
        size_t pos = find_in_row(i, j);
        if (pos != static_cast<size_t>(-1)) return data_[pos];
        return T{};
    }

    // Set element at position (i,j), updating CSR structure as needed
    void set(size_t i, size_t j, const T &value) {
        if (i >= nrows_ || j >= ncols_) throw invalid_index();
        size_t start = indptr_[i];
        size_t end = indptr_[i + 1];
        // Try find existing
        // Binary search in [start, end)
        size_t l = start, r = end;
        while (l < r) {
            size_t mid = l + (r - l) / 2;
            if (indices_[mid] == j) {
                data_[mid] = value;
                return;
            }
            if (indices_[mid] < j) l = mid + 1; else r = mid;
        }
        // Insert new at position l (first >= j) to keep sorted
        size_t pos = l;
        indices_.insert(indices_.begin() + pos, j);
        data_.insert(data_.begin() + pos, value);
        for (size_t row = i + 1; row < indptr_.size(); ++row) {
            ++indptr_[row];
        }
    }

    // Return the row pointer array
    const std::vector<size_t> &getIndptr() const { return indptr_; }

    // Return the column indices array
    const std::vector<size_t> &getIndices() const { return indices_; }

    // Return the data values array
    const std::vector<T> &getData() const { return data_; }

    // Convert CSR format to dense matrix representation
    std::vector<std::vector<T>> getMatrix() const {
        std::vector<std::vector<T>> dense(nrows_, std::vector<T>(ncols_, T{}));
        for (size_t i = 0; i < nrows_; ++i) {
            for (size_t k = indptr_[i]; k < indptr_[i + 1]; ++k) {
                size_t j = indices_[k];
                dense[i][j] = data_[k];
            }
        }
        return dense;
    }

    // Matrix-vector multiplication
    std::vector<T> operator*(const std::vector<T> &vec) const {
        if (vec.size() != ncols_) throw size_mismatch();
        std::vector<T> res(nrows_, T{});
        for (size_t i = 0; i < nrows_; ++i) {
            T acc{};
            for (size_t k = indptr_[i]; k < indptr_[i + 1]; ++k) {
                acc = acc + data_[k] * vec[indices_[k]];
            }
            res[i] = acc;
        }
        return res;
    }

    // Extract submatrix containing rows [l,r)
    CSRMatrix getRowSlice(size_t l, size_t r) const {
        if (l > r || r > nrows_) throw invalid_index();
        CSRMatrix sub(r - l, ncols_);
        size_t start = indptr_[l];
        size_t end = indptr_[r];
        size_t cnt = end - start;
        sub.indices_.assign(indices_.begin() + start, indices_.begin() + end);
        sub.data_.assign(data_.begin() + start, data_.begin() + end);
        sub.indptr_.resize(sub.nrows_ + 1);
        sub.indptr_[0] = 0;
        for (size_t i = 0; i < sub.nrows_; ++i) {
            sub.indptr_[i + 1] = indptr_[l + i + 1] - start;
        }
        return sub;
    }
};

}

#endif // CSR_MATRIX_HPP
