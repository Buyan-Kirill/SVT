#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

double f(double x)
{
    return -30 * pow(x, 4);
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

void fill_right_part(std::vector<double> &f_vec)
{
    int N = f_vec.size();
    double h = 1. / (N + 1);
    for (int i = 1; i <= N; ++i) {
        f_vec[i - 1] = f(h * i) * h * h;
    }
    f_vec[N - 1] += 1.;
}

int main()
{
    std::ofstream myfile;
    myfile.open ("output.txt");
    int n = 5;
    double a = 0., b = 1.;
    for (int i = 0; i < 9; i++) {
        std::vector<double> A(n * n), f_vec(n), y(n + 2);

        fill_matrix(A, n);
        fill_right_part(f_vec);
        
        progonka(A, f_vec, y.end() - 2);
        y[0] = a;
        y[n + 1] = b;
  
        for (auto i: y) {
            myfile << i << " ";
        }
        myfile << '\n';
        n *= 3;
    }
    myfile.close();
}
