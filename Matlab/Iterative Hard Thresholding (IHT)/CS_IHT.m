
clc
clear all

% Iterative Hard Thresholding (IHT) implementation

% generate sparse signal x of length n = 1000 with 100 non-zero coefficients
n = 1000;
s = 100; % sparsity
sigma = 1;
x = zeros(n, 1);
nz = randperm(n, s); % randomly choose 100 indices
x(nz) = sigma*randn(s, 1); % set non-zero coefficients to Gaussian noise

% generate observation matrix A
m = 500; % number of measurements
A = randn(m, n);

% generate noisy observation vector b
sigma = 0.1; % noise standard deviation
b = A*x + sigma*randn(m, 1);

% set IHT parameters
lambda = -0.999992; % regularization parameter
tol = 1e-6; % convergence tolerance
maxit = 1000; % maximum number of iterations

% run IHT algorithm
xhat = iht(A, b, lambda, tol, maxit);

% Compute metrics to evaluate the quality of the recovered signal
mse = immse(xhat, x);
psnr = psnr(xhat, x);
ssimval = ssim(x, xhat);
sqs = round(10*log10(ssimval^2)); % subjective quality scale

% display results
fprintf('MSE: %f\n', mse);
fprintf('PSNR: %f dB\n', psnr);
fprintf('SSIM: %f\n', ssimval);
fprintf('SQS = %d\n', sqs);
% plot x and xhat
figure;
plot(x, 'b'); hold on;
plot(xhat, 'r--'); hold off;
xlabel('Index'); ylabel('Value');
legend('True signal', 'Recovered signal');
title('IHT')

% IHT function implementation
function x = iht(A, b, lambda, tol, maxit)
    % initialize variables
    x = zeros(size(A, 2), 1);
    k = 0;
    
    % main iteration loop
    while k < maxit
        k = k + 1;
        xold = x;
        x = shrink(A'*b + lambda*x, lambda);
        if norm(x - xold)/norm(x) < tol
            break;
        end
    end
end

% shrinkage operator implementation
function y = shrink(x, lambda)
    y = x.*(abs(x) >= lambda);
end
