#pragma once

#include <vector>

/** @brief One sparse matrix entry before CSR compression. */
struct SparseTriplet { int row = -1; int column = -1; double value = 0.0; };

/** @brief Deterministic compressed-sparse-row matrix. */
class CsrMatrix {
public:
    static CsrMatrix fromTriplets(int rows, int columns, std::vector<SparseTriplet> triplets, double numerical_epsilon = 1e-18);
    int rows() const;
    int columns() const;
    void multiply(const std::vector<double>& x, std::vector<double>& result) const;
    std::vector<double> diagonal() const;
private:
    int rows_ = 0, columns_ = 0;
    std::vector<int> row_ptr_;
    std::vector<int> column_indices_;
    std::vector<double> values_;
};
