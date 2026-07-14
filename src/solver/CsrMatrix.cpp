#include "solver/CsrMatrix.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

CsrMatrix CsrMatrix::fromTriplets(int rows, int columns, std::vector<SparseTriplet> triplets, double eps) {
    if (rows < 0 || columns < 0) throw std::invalid_argument("CsrMatrix::fromTriplets: negative dimensions");
    if (!std::isfinite(eps) || eps < 0.0) throw std::invalid_argument("CsrMatrix::fromTriplets: invalid epsilon");
    for (const auto& t : triplets) {
        if (t.row < 0 || t.row >= rows || t.column < 0 || t.column >= columns) throw std::out_of_range("CsrMatrix::fromTriplets: triplet index out of range");
        if (!std::isfinite(t.value)) throw std::invalid_argument("CsrMatrix::fromTriplets: non-finite value");
    }
    std::sort(triplets.begin(), triplets.end(), [](const auto& a, const auto& b){ return a.row != b.row ? a.row < b.row : a.column < b.column; });
    CsrMatrix m; m.rows_ = rows; m.columns_ = columns; m.row_ptr_.assign(static_cast<size_t>(rows)+1, 0);
    for (size_t i=0;i<triplets.size();) {
        int r=triplets[i].row, c=triplets[i].column; double v=0.0;
        while (i<triplets.size() && triplets[i].row==r && triplets[i].column==c) v += triplets[i++].value;
        if (std::abs(v) > eps) { ++m.row_ptr_[static_cast<size_t>(r)+1]; m.column_indices_.push_back(c); m.values_.push_back(v); }
    }
    for (int r=0;r<rows;++r) m.row_ptr_[static_cast<size_t>(r)+1] += m.row_ptr_[static_cast<size_t>(r)];
    return m;
}
int CsrMatrix::rows() const { return rows_; }
int CsrMatrix::columns() const { return columns_; }
void CsrMatrix::multiply(const std::vector<double>& x, std::vector<double>& result) const {
    if (static_cast<int>(x.size()) != columns_) throw std::invalid_argument("CsrMatrix::multiply: dimension mismatch");
    result.assign(static_cast<size_t>(rows_), 0.0);
    for (int r=0;r<rows_;++r) for (int k=row_ptr_[r]; k<row_ptr_[r+1]; ++k) result[static_cast<size_t>(r)] += values_[static_cast<size_t>(k)] * x[static_cast<size_t>(column_indices_[static_cast<size_t>(k)])];
}
std::vector<double> CsrMatrix::diagonal() const {
    std::vector<double> d(static_cast<size_t>(rows_), 0.0);
    for (int r=0;r<rows_;++r) for (int k=row_ptr_[r]; k<row_ptr_[r+1]; ++k) if (column_indices_[static_cast<size_t>(k)] == r) d[static_cast<size_t>(r)] += values_[static_cast<size_t>(k)];
    return d;
}
