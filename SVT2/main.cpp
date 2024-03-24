#include "inmost.h"
#include <stdio.h>

using namespace INMOST;


void fill_matrix(Sparse::Matrix &A)
{
    const int N = A.Size();
    int M = std::sqrt(N);
    for (int i = 0; i < N; ++i) {
        int y = i / M;
        int x = i % M;
        if (y + 1 < M) {
            A[i][(y + 1) * M + x] = -1.;
        }
        if (y > 0) {
            A[i][(y - 1) * M + x] = -1.;
        }
        if (x + 1 < M) {
            A[i][y * M + x + 1] = -1.;
        }
        if (x > 0) {
            A[i][y * M + x - 1] = -1.;
        }
        A[i][i] = 4.;
    }
}

double u(double x, double y)
{
    return std::sin(5 * x) * std::cos(y);
}

double f(double x, double y)
{
    return 26 * std::sin(5 * x) * std::cos(y);
}

void fill_right_part(Sparse::Vector &f_vec)
{
    int N = f_vec.Size();
    int M = std::sqrt(N);
    double h = 1. / (M + 1);
    for (int i = 0; i < N; ++i) {
        int y = i / M;
        int x = i % M;
        f_vec[i] = f((x + 1) * h, (y + 1) * h) * h * h;
        if (y == M - 1) {
            f_vec[i] += u((x + 1) * h, 1.);
        }
        if (y == 0) {
            f_vec[i] += u((x + 1) * h, 0.);
        }
        if (x == M - 1) {
            f_vec[i] += u(1., (y + 1) * h);
        }
        if (x == 0) {
            f_vec[i] += u(0., (y + 1) * h);
        }
    }
}

double C_norm (double (*u)(double x, double y), double a, double b, Sparse::Vector &sol)
{
    double c_norm = 0.;
    int N = sol.Size();
    int M = sqrt(N);
    double h = (b - a) / (M + 1);
    for (int i = 0; i < N; i++) {
        int y = i / M;
        int x = i % M;
        c_norm = std::max(c_norm, fabs(u((x + 1) * h, (y + 1) * h) - sol[i]));
    }
    return c_norm;
}

double L_norm (double (*u)(double x, double y), double a, double b, Sparse::Vector &sol)
{
    double l_norm = 0.;
    int N = sol.Size();
    int M = sqrt(N);
    double h = (b - a) / (M + 1);
    for (int i = 0; i < N; i++) {
        int y = i / M;
        int x = i % M;
        l_norm += (u((x + 1) * h, (y + 1) * h) - sol[i]) * (u((x + 1) * h, (y + 1) * h) - sol[i]) * h;
    }
    return sqrt(l_norm);
}

void printer (const std::vector<double> &a)
{
    std::cout << "[";
    for (int i = 0; i < a.size(); i++) {
        if(i > 0)
            std::cout << ", ";
        std::cout << a[i];
    }
    std::cout << "]" << std::endl;
}

int main(int argc, char *argv[])
{
    std::vector<double> c_errors;
    std::vector<double> l_errors;
    std::vector<double> iter_times;
    std::vector<double> times;

    for (int N = 8; N <= 128; N *= 2) {
        // Create sparse matrix, RHS vector and solution vector
        Sparse::Matrix A;
        Sparse::Vector b;
        Sparse::Vector sol;
        // Set their size
        A.SetInterval(0, N * N);
        b.SetInterval(0, N * N);
        sol.SetInterval(0, N * N);

        // std::cout << "Calculating matrix\n";
        fill_matrix(A);
        // std::cout << "Calculating right side\n";
        fill_right_part(b);

        // printer(A);
        // printer(b);
        

        // Get solver
        // All inner INMOST solvers are BiCGStab
        // with different preconditioners, let's use ILU2
        // std::cout << "Solving system\n";
        Solver S(Solver::INNER_DDPQILUC);
        S.SetParameter("absolute_tolerance", "1e-12");
        S.SetParameter("relative_tolerance", "1e-10");

        // Set matrix in the solver;
        // this also computes preconditioner
        S.SetMatrix(A);

        // Solve
        bool solved = S.Solve(b, sol);
        // std::cout << "prec.time: " << S.PreconditionerTime() << std::endl;
        iter_times.emplace_back(S.IterationsTime());
        times.emplace_back(S.Iterations() * S.IterationsTime());
        if(!solved){
            std::cout << "Linear solver failure!" << std::endl;
            std::cout << "Reason: " << S.ReturnReason() << std::endl;
        }

        l_errors.emplace_back(L_norm(u, 0., 1., sol));
        c_errors.emplace_back(C_norm(u, 0., 1., sol));
        std::cout << N << std::endl;
    }
    
    printer(c_errors);
    printer(l_errors);
    printer(times);
    printer(iter_times);

	return 0;
}
