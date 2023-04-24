% shakya_jayakody_hw4

clc
clear all

% Orthogonal Matching Pursuit (OMP) implementation

% generate sparse signal x of length n = 1000 with 100 non-zero coefficients
n = 1000;
x = zeros(n, 1);
nz = randsample(n, 100); % randomly choose 100 indices
x(nz) = randn(100, 1); % set non-zero coefficients to Gaussian noise

% generate observation matrix A
m = 500; % number of measurements
A = randn(m, n);

% generate noisy observation vector b
sigma = 0.1; % noise standard deviation
b = A*x + sigma*randn(m, 1);

% set OMP parameters
k = 100; % sparsity level
tol = 1e-6; % convergence tolerance

% run OMP algorithm
xhat = omp(A, b, k, tol);

% compute MSE, PSNR, SSIM between x and xhat
mse = norm(x - xhat)^2/n;
psnr = 20*log10(max(abs(x))/sqrt(mse));
ssimval = ssim(x, xhat);

% display results
fprintf('MSE: %.4f\n', mse);
fprintf('PSNR: %.4f dB\n', psnr);
fprintf('SSIM: %.4f\n', ssimval);

% plot x and xhat
figure;
plot(x, 'b'); hold on;
plot(xhat, 'r--'); hold off;
xlabel('Index'); ylabel('Value');
legend('True signal', 'Recovered signal');
title('OMP')

% OMP function implementation
function x = omp(A, b, k, tol)
    % initialize variables
    x = zeros(size(A, 2), 1);
    r = b;
    T = [];
    
    % main iteration loop
    while length(T) < k
        % compute correlation between residual and columns of A
        corrs = abs(A'*r);
        % add index of highest correlation to support set
        T = [T; find(corrs == max(corrs), 1)];
        % solve least squares problem for support set
        x(T) = A(:, T)\b;
        % update residual
        r = b - A*x;
        if norm(r) < tol
            break;
        end
    end
end
