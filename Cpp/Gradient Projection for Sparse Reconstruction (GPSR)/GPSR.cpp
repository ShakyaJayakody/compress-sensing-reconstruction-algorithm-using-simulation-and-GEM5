#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> gpsr(vector<vector<double>> A, vector<double> b, double lambda, int max_iters)
{
    int m = A.size();
    int n = A[0].size();
    vector<double> x(n, 0);
    vector<double> r(m);
    vector<double> w(n);
    double tau = 2;
    double beta = 1;
    
    for(int i = 0; i < max_iters; i++) {
        // Compute gradient
        for(int j = 0; j < n; j++) {
            double dot = 0;
            for(int k = 0; k < m; k++) {
                dot += A[k][j] * r[k];
            }
            w[j] = x[j] - dot;
        }
        // Compute thresholding operator
        double alpha = lambda / beta;
        for(int j = 0; j < n; j++) {
            if(abs(w[j]) <= alpha) {
                x[j] = 0;
            } else if(w[j] > alpha) {
                x[j] = w[j] - alpha;
            } else {
                x[j] = w[j] + alpha;
            }
        }
        // Compute projection
        double norm_x = 0;
        for(int j = 0; j < n; j++) {
            norm_x += x[j] * x[j];
        }
        norm_x = sqrt(norm_x);
        if(norm_x > 1) {
            for(int j = 0; j < n; j++) {
                x[j] /= norm_x;
            }
        }
        // Compute residual
        for(int k = 0; k < m; k++) {
            double dot = 0;
            for(int j = 0; j < n; j++) {
                dot += A[k][j] * x[j];
            }
            r[k] = b[k] - dot;
        }
        // Update beta and tau
        beta = norm_x;
        tau *= beta / norm_x;
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
    vector<double> x = gpsr(A, b, lambda, max_iters);
    // Use x for sparse signal reconstruction
    // ...
    return 0;
}

