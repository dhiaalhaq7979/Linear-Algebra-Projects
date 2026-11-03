#include "Matrix.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

static const double EPS = 1e-12;

Matrix::Matrix(int rows, int cols, double fill) : r(rows), c(cols), a(rows, std::vector<double>(cols, fill)) {}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> values)
    : r(static_cast<int>(values.size())), c(static_cast<int>(values.begin()->size())) {
    for (const auto& row : values) {
        if (static_cast<int>(row.size()) != c) throw std::invalid_argument("all rows must have the same length");
        a.emplace_back(row);
    }
}

Matrix Matrix::identity(int n) {
    Matrix m(n, n);
    for (int i = 0; i < n; i++) m.a[i][i] = 1;
    return m;
}

Matrix Matrix::operator+(const Matrix& o) const {
    if (r != o.r || c != o.c) throw std::invalid_argument("addition needs matrices of the same size");
    Matrix m(r, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) m.a[i][j] = a[i][j] + o.a[i][j];
    return m;
}

Matrix Matrix::operator-(const Matrix& o) const { return *this + (o * -1.0); }

Matrix Matrix::operator*(const Matrix& o) const {
    if (c != o.r) throw std::invalid_argument("multiplication needs cols(A) == rows(B)");
    Matrix m(r, o.c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < o.c; j++)
            for (int k = 0; k < c; k++) m.a[i][j] += a[i][k] * o.a[k][j];
    return m;
}

Matrix Matrix::operator*(double k) const {
    Matrix m(*this);
    for (auto& row : m.a)
        for (double& x : row) x *= k;
    return m;
}

bool Matrix::approxEquals(const Matrix& o, double eps) const {
    if (r != o.r || c != o.c) return false;
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            if (std::fabs(a[i][j] - o.a[i][j]) > eps) return false;
    return true;
}

Matrix Matrix::transpose() const {
    Matrix t(c, r);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) t.a[j][i] = a[i][j];
    return t;
}

double Matrix::determinantCofactor() const {
    if (r != c) throw std::invalid_argument("determinant needs a square matrix");
    if (r == 1) return a[0][0];
    if (r == 2) return a[0][0] * a[1][1] - a[0][1] * a[1][0];
    double det = 0;
    for (int col = 0; col < c; col++) {
        Matrix minor(r - 1, c - 1);
        for (int i = 1; i < r; i++)
            for (int j = 0, mj = 0; j < c; j++)
                if (j != col) minor.a[i - 1][mj++] = a[i][j];
        det += (col % 2 == 0 ? 1 : -1) * a[0][col] * minor.determinantCofactor();
    }
    return det;
}

double Matrix::determinant() const {
    if (r != c) throw std::invalid_argument("determinant needs a square matrix");
    Matrix m(*this);
    double det = 1;
    for (int col = 0; col < r; col++) {
        int pivot = col;
        for (int i = col + 1; i < r; i++)
            if (std::fabs(m.a[i][col]) > std::fabs(m.a[pivot][col])) pivot = i;
        if (std::fabs(m.a[pivot][col]) < EPS) return 0;
        if (pivot != col) { std::swap(m.a[pivot], m.a[col]); det = -det; }   // a row swap flips the sign
        det *= m.a[col][col];
        for (int i = col + 1; i < r; i++) {
            double f = m.a[i][col] / m.a[col][col];
            for (int j = col; j < c; j++) m.a[i][j] -= f * m.a[col][j];
        }
    }
    return det;
}

Matrix Matrix::inverse() const {
    if (r != c) throw std::invalid_argument("inverse needs a square matrix");
    Matrix m(*this), inv = identity(r);
    for (int col = 0; col < r; col++) {
        int pivot = col;
        for (int i = col + 1; i < r; i++)
            if (std::fabs(m.a[i][col]) > std::fabs(m.a[pivot][col])) pivot = i;
        if (std::fabs(m.a[pivot][col]) < EPS) throw std::runtime_error("matrix is singular (no inverse)");
        std::swap(m.a[pivot], m.a[col]);
        std::swap(inv.a[pivot], inv.a[col]);
        double p = m.a[col][col];
        for (int j = 0; j < r; j++) { m.a[col][j] /= p; inv.a[col][j] /= p; }
        for (int i = 0; i < r; i++) {
            if (i == col) continue;
            double f = m.a[i][col];
            for (int j = 0; j < r; j++) { m.a[i][j] -= f * m.a[col][j]; inv.a[i][j] -= f * inv.a[col][j]; }
        }
    }
    return inv;
}

int Matrix::rank() const {
    Matrix m(*this);
    int rank = 0;
    for (int col = 0; col < c && rank < r; col++) {
        int pivot = rank;
        for (int i = rank + 1; i < r; i++)
            if (std::fabs(m.a[i][col]) > std::fabs(m.a[pivot][col])) pivot = i;
        if (std::fabs(m.a[pivot][col]) < 1e-9) continue;
        std::swap(m.a[pivot], m.a[rank]);
        for (int i = rank + 1; i < r; i++) {
            double f = m.a[i][col] / m.a[rank][col];
            for (int j = col; j < c; j++) m.a[i][j] -= f * m.a[rank][j];
        }
        rank++;
    }
    return rank;
}

void Matrix::print(const std::string& title, int precision) const {
    std::cout << title << " (" << r << "x" << c << ")\n";
    for (const auto& row : a) {
        std::cout << "  [";
        for (int j = 0; j < c; j++) {
            double v = std::fabs(row[j]) < 1e-10 ? 0.0 : row[j];
            std::cout << std::setw(precision + 6) << std::fixed << std::setprecision(precision) << v;
        }
        std::cout << " ]\n";
    }
}

std::vector<double> solveLinearSystem(const Matrix& A, const std::vector<double>& b) {
    int n = A.rows();
    if (A.cols() != n || static_cast<int>(b.size()) != n) throw std::invalid_argument("A must be n x n and b of size n");
    // augmented matrix [A | b]
    std::vector<std::vector<double>> m(n, std::vector<double>(n + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) m[i][j] = A.at(i, j);
        m[i][n] = b[i];
    }
    for (int col = 0; col < n; col++) {                 // forward elimination
        int pivot = col;
        for (int i = col + 1; i < n; i++)
            if (std::fabs(m[i][col]) > std::fabs(m[pivot][col])) pivot = i;
        if (std::fabs(m[pivot][col]) < EPS) throw std::runtime_error("system has no unique solution");
        std::swap(m[pivot], m[col]);
        for (int i = col + 1; i < n; i++) {
            double f = m[i][col] / m[col][col];
            for (int j = col; j <= n; j++) m[i][j] -= f * m[col][j];
        }
    }
    std::vector<double> x(n);                           // back substitution
    for (int i = n - 1; i >= 0; i--) {
        double s = m[i][n];
        for (int j = i + 1; j < n; j++) s -= m[i][j] * x[j];
        x[i] = s / m[i][i];
    }
    return x;
}
