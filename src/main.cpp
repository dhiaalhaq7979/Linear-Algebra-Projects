// Linear Algebra demo - every operation of the Matrix class, with checks
#include <cmath>
#include <iomanip>
#include <iostream>

#include "Matrix.h"

int main() {
    std::cout << std::fixed << std::setprecision(3);
    Matrix A = {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};
    Matrix B = {{1, 0, 2}, {0, 3, -1}, {4, 1, 0}};

    std::cout << "===== Basic operations =====\n";
    A.print("A");
    B.print("B");
    (A + B).print("A + B");
    (A - B).print("A - B");
    (A * B).print("A x B");
    (B * A).print("B x A");
    std::cout << "A x B == B x A ? " << ((A * B).approxEquals(B * A) ? "yes" : "no (matrix multiplication is not commutative)") << "\n\n";
    A.transpose().print("transpose(A)");
    std::cout << "(A x B)^T == B^T x A^T ? " << ((A * B).transpose().approxEquals(B.transpose() * A.transpose()) ? "yes" : "no") << "\n\n";

    std::cout << "===== Determinant =====\n";
    Matrix C2 = {{4, 7}, {2, 6}};
    std::cout << "det([[4,7],[2,6]])           = " << C2.determinant() << '\n';
    std::cout << "det(A) by cofactor expansion = " << A.determinantCofactor() << '\n';
    std::cout << "det(A) by Gaussian elimination = " << A.determinant() << '\n';
    Matrix D4 = {{1, 2, 0, 3}, {0, 1, 4, 1}, {2, 0, 1, 0}, {1, 1, 1, 1}};
    std::cout << "4x4 matrix: cofactor = " << D4.determinantCofactor() << ", elimination = " << D4.determinant() << '\n';
    std::cout << "det(A x B) = " << (A * B).determinant() << ",  det(A) x det(B) = " << A.determinant() * B.determinant() << "\n\n";

    std::cout << "===== Inverse =====\n";
    Matrix Ainv = A.inverse();
    Ainv.print("inverse(A)");
    (A * Ainv).print("A x inverse(A)");
    std::cout << "equals identity? " << ((A * Ainv).approxEquals(Matrix::identity(3), 1e-9) ? "yes" : "no") << '\n';
    Matrix S = {{1, 2, 3}, {2, 4, 6}, {1, 0, 1}};
    std::cout << "Singular matrix (row 2 = 2 x row 1): det = " << S.determinant() << ", rank = " << S.rank() << '\n';
    try { S.inverse(); } catch (const std::exception& e) { std::cout << "inverse -> error: " << e.what() << "\n\n"; }

    std::cout << "===== Solving A x = b (Gaussian elimination) =====\n";
    std::cout << "   2x +  y -  z =   8\n  -3x -  y + 2z = -11\n  -2x +  y + 2z =  -3\n";
    std::vector<double> b = {8, -11, -3};
    std::vector<double> x = solveLinearSystem(A, b);
    std::cout << "Solution: x = " << x[0] << ", y = " << x[1] << ", z = " << x[2] << '\n';
    std::cout << "Check A x = ";
    for (int i = 0; i < 3; i++) {
        double s = 0;
        for (int j = 0; j < 3; j++) s += A.at(i, j) * x[j];
        std::cout << s << (i < 2 ? ", " : "\n");
    }
    Matrix E = {{0, 2, 1, -1}, {1, 1, 0, 2}, {3, 0, -1, 1}, {2, -1, 3, 0}};
    std::vector<double> x4 = solveLinearSystem(E, {3, 7, 4, 8});
    std::cout << "4x4 system (needs a row swap because E[0][0] = 0): ";
    for (double v : x4) std::cout << v << ' ';
    double worst = 0;
    std::vector<double> rhs = {3, 7, 4, 8};
    for (int i = 0; i < 4; i++) {
        double s = 0;
        for (int j = 0; j < 4; j++) s += E.at(i, j) * x4[j];
        worst = std::fmax(worst, std::fabs(s - rhs[i]));
    }
    std::cout << "\nlargest residual |E x - b| = " << std::scientific << std::setprecision(1) << worst << '\n';
    return 0;
}
