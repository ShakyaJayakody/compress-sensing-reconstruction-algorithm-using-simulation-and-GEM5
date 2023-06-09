#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> iht(vector<vector<double>> A, vector<double> b, double lambda, int max_iters)
{
    int m = A.size();
    int n = A[0].size();
    vector<double> x(n, 0);
    
    for(int i = 0; i < max_iters; i++) {
            // Compute residual
        vector<double> r(m, 0);
        for(int k = 0; k < m; k++) {
            for(int j = 0; j < n; j++) {
                r[k] += A[k][j] * x[j];
            }
            r[k] -= b[k];
        }
        
        // Compute gradient
        vector<double> grad(n, 0);
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < m; k++) {
                grad[j] += A[k][j] * r[k];
            }
        }
        
        // Apply hard thresholding
        for(int j = 0; j < n; j++) {
            if(abs(x[j] + grad[j] / lambda) <= 1 / lambda) {
                x[j] = 0;
            } else {
                x[j] = x[j] + grad[j] / lambda;
            }
        }
    }
    
    return x;
}

int main()
{
    int n = 100;
    int m = 50;
    double lambda = 0.1;
    int max_iters = 1000;
    vector<vector<double>> A(m, vector<double>(n));
    vector<double> b(m);
    // Initialize A and b
    // ...
    vector<double> x = iht(A, b, lambda, max_iters);
    // Use x for sparse signal reconstruction
    // ...
    return 0;
}
