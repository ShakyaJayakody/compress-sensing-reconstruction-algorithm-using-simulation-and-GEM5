% shakya_jayakody_hw4

clc
clear all

% Generate sparse signal x
n = 1000; % signal length
s = 100; % sparsity
sigma = 1; % standard deviation of non-zero coefficients
x = zeros(n, 1);
nz_idx = randperm(n, s); % indices of non-zero coefficients
x(nz_idx) = sigma * randn(s, 1); % non-zero coefficients from Gaussian distribution

% Add noise to signal
noise_sigma = 0.5; % standard deviation of noise
noise = noise_sigma * randn(size(x));
y = x + noise;

% Set parameters for IST
lambda = 0.1; % regularization parameter
max_iter = 100; % maximum number of iterations

% Define soft thresholding function
soft_thresh = @(x, T) sign(x) .* max(abs(x) - T, 0);

% Run IST algorithm
x_hat = y; % initial estimate
% plot x and xhat
figure;

plot(x, 'b'); hold on;
plot(x_hat, 'r--'); hold off;
xlabel('Index'); ylabel('Value');
legend('True signal', 'Recovered signal');
title('IHT')

for i = 1:max_iter
    x_hat_old = x_hat;
    x_hat = soft_thresh(x_hat + lambda * (y - x_hat), lambda);
    % Check for convergence
    if norm(x_hat - x_hat_old) < 1e-6
        break;
    end
end

% Compute metrics to evaluate the quality of the denoised signal
mse = immse(x_hat, x);
psnr = psnr(x_hat, x);
ssim_val = ssim(x_hat, x);
sqs = round(10*log10(ssim_val^2)); % subjective quality scale
fprintf('MSE = %f\n', mse);
fprintf('PSNR = %f dB\n', psnr);
fprintf('SSIM = %f\n', ssim_val);
fprintf('SQS = %d\n', sqs);

