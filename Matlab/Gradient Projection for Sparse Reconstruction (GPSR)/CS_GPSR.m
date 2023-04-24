% shakya_jayakody_hw4

clc
clear all
% Gradient Projection for Sparse Reconstruction (GPSR) implementation

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

% set GPSR parameters
lambda = 0.1; % regularization parameter
tol = 1e-6; % convergence tolerance
maxit = 1000; % maximum number of iterations

% run GPSR algorithm
xhat = gpsr(A, b, lambda, tol, maxit);

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
title('GPSR')

% GPSR function implementation
function x = gpsr(A, b, lambda, tol, maxit)
    % initialize variables
    x = zeros(size(A, 2), 1);
    r = b;
    w = x;
    v = zeros(size(A, 1), 1);
    k = 0;
    mu = norm(A'*r, 'inf');
    rho = 1.1;
    
    % main iteration loop
    while k < maxit && norm(r)/norm(b) > tol
        k = k + 1;
        v = A*w - b;
        u = A'*v;
        xold = x;
        x = shrink(w - (1/mu)*u, lambda/mu);
        w = x + ((k-2)/(k+1))*(x - xold);
        mu = rho*mu;
    end
end

% shrinkage operator implementation
function y = shrink(x, lambda)
    y = sign(x).*max(abs(x) - lambda, 0);
end
