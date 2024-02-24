#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

double u(double x)
{
    return -30 * pow(x, 4);
}

void printer(std::vector<double> y)
{
    for (auto i: y) {
        std::cout << i << ", ";
    }
    std::cout << std::endl << std::endl;
}

void mat_printer(std::vector<double> A, int N)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << A[i * N + j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void progonka(std::vector<double> A, std::vector<double> f, std::vector<double>::iterator it)
{
    const int N = f.size();

    for (int i = 0; i < N - 1; i++) {
        A[(i + 1) * N + i + 1] -= A[i * N + i + 1] * A[(i + 1) * N + i] / A[i * N + i];
        f[i + 1] -= f[i] * A[(i + 1) * N + i] / A[i * N + i];
        A[(i + 1) * N + i] = 0;
    }

    *it = f[N - 1] / A[N * N - 1];
    for (int i = N - 2; i >= 0; i--) {
        auto last = *it;
        --it;
        *it = (f[i] - last * A[i * N + i + 1]) / A[i * N + i];
    }
}

void mat_vec(const std::vector<double> &A, const std::vector<double> &x, std::vector<double> &y)
{
    for (size_t i = 0; i < x.size(); i++) {
        y[i] = 0;
        for (size_t j = 0; j < x.size(); j++) {
            y[i] += A[i * x.size() + j] * x[j];
        }
    }
}

void fill_matrix(std::vector<double> &A, const int &N)
{
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i * N + j] = 0;
            if (abs(i - j) <= 1) {
                A[i * N + j] = -1;
            } 
            if (i == j) {
                A[i * N + j] = 2;
            }
        }
    }
}

void fill_right_part(std::vector<double> &f)
{
    int N = f.size();
    double h = 1. / (N + 1);
    for (int i = 1; i <= N; ++i) {
        f[i - 1] = u(h * i) * h * h;
    }
    f[N - 1] += 1.;
}

int main()
{
    std::ofstream myfile;
    myfile.open ("output.txt");
    int n = 5;
    for (int i = 0; i < 9; i++) {
        std::vector<double> A(n * n), f(n), x(n + 2);

        fill_matrix(A, n);
        fill_right_part(f);
        
        progonka(A, f, x.end() - 2);
        x[0] = 0;
        x[n + 1] = 1;
  
        for (auto i: x) {
            myfile << i << " ";
        }
        myfile << '\n';
        n *= 3;
    }
    myfile.close();
}
