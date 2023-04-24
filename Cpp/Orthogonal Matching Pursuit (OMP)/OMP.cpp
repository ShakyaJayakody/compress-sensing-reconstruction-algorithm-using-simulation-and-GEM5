#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> omp(vector<vector<double>> A, vector<double> b, int k)
{
    int m = A.size();
    int n = A[0].size();
    vector<double> x(n, 0);
    vector<double> r = b;
    vector<int> omega(k);
    int omega_size = 0;
    
    while(omega_size < k) {
        // Compute correlations
        vector<double> corr(n, 0);
        for(int j = 0; j < n; j++) {
            if(x[j] != 0) {
                continue;
            }
            double dot = 0;
            for(int k = 0; k < m; k++) {
                dot += A[k][j] * r[k];
            }
            corr[j] = abs(dot);
        }
        // Select new index
        int j = distance(corr.begin(), max_element(corr.begin(), corr.end()));
        omega[omega_size++] = j;
        // Solve least-squares problem
        vector<vector<double>> At(omega_size, vector<double>(m));
        for(int i = 0; i < omega_size; i++) {
            for(int k = 0; k < m; k++) {
                At[i][k] = A[k][omega[i]];
            }
        }
        vector<double> y(omega_size);
        for(int i = 0; i < omega_size; i++) {
            y[i] = r[i];
        }
        vector<double> xk(omega_size);
        for(int i = 0; i < omega_size; i++) {
            xk[i] = x[omega[i]];
        }
        vector<double> zk = lsq(At, y);
        // Update x and residual
        for(int i = 0; i < omega_size; i++) {
            x[omega[i]] = zk[i];
        }
        for(int k = 0; k < m; k++) {
            double dot = 0;
            for(int j = 0; j < n; j++) {
                dot += A[k][j] * x[j];
            }
            r[k] = b[k] - dot;
        }
    }
    
    return x;
}

vector<double> lsq(vector<vector<double>> A, vector<double> b)
{
    int m = A.size();
    int n = A[0].size();
    vector<vector<double>> Q(m, vector<double>(n));
    vector<vector<double>> R(n, vector<double>(n));
    // QR factorization
    for(int j = 0; j < n; j++) {
        for(int i = 0; i < m; i++) {
            Q[i][j] = A[i][j];
        }
        for(int k = 0; k < j; k++) {
            double dot = 0;
            for(int i = 0; i < m; i++) {
                dot += Q[i][j] * Q[i][k];
            }
            R[k][j] = dot;
            for(int i = 0; i < m; i++) {
                Q[i][j] -= R[k][j] * Q[i][
				
			}

int main()
{
  
    // ...
    return 0;
}
