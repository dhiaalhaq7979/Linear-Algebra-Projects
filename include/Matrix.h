#ifndef MATRIX_H
#define MATRIX_H

#include <initializer_list>
#include <string>
#include <vector>

// A small dense matrix class for linear algebra practice.
class Matrix {
private:
    int r, c;
    std::vector<std::vector<double>> a;

public:
    Matrix(int rows, int cols, double fill = 0.0);
    Matrix(std::initializer_list<std::initializer_list<double>> values);
    static Matrix identity(int n);

    int rows() const { return r; }
    int cols() const { return c; }
    double& at(int i, int j) { return a[i][j]; }
    double at(int i, int j) const { return a[i][j]; }

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double k) const;
    bool approxEquals(const Matrix& other, double eps = 1e-9) const;

    Matrix transpose() const;
    double determinantCofactor() const;   // Laplace expansion (fine for small n)
    double determinant() const;           // Gaussian elimination, O(n^3)
    Matrix inverse() const;               // Gauss-Jordan, throws if singular
    int rank() const;

    void print(const std::string& title, int precision = 3) const;
};

// Solve A x = b with Gaussian elimination and partial pivoting.
std::vector<double> solveLinearSystem(const Matrix& A, const std::vector<double>& b);

#endif
