%%
close all;
clc;
clear;


%% Signal model
% Define parameters
M = 10; % Number of antennas
N = 20; % Number of samples
d = 2; % Number of sources
Delta = 0.5; % Antenna spacing in wavelengths
theta = [10, -15]'; % Directions of the sources in degrees
f = [0.3, 0.1]'; % Normalized frequencies of the sources
SNR = 100; % Signal-to-noise ratio per source in dB

% Generate data
[X, ~, ~] = gendata(M, N, Delta, theta, f, SNR);
%%
% Compute and plot singular values
singular_values = svd(X);
figure;
stem(singular_values, 'filled');
title('Singular Values of the Received Signal Matrix X');
xlabel('Index');
ylabel('Singular Value');
% saveas(gcf,'singular.png')

% Scenario 1: Double the number of samples
N2 = 2 * N;
[X2, ~, ~] = gendata(M, N2, Delta, theta, f, SNR);
singular_values2 = svd(X2);
figure;
stem(singular_values2, 'filled');
title('Singular Values of X with Doubled Number of Samples');
xlabel('Index');
ylabel('Singular Value');
% saveas(gcf,'singulardoubleN.png')

% Scenario 2: Double the number of antennas
M2 = 2 * M;
[X3, ~, ~] = gendata(M2, N, Delta, theta, f, SNR);
singular_values3 = svd(X3);
figure;
stem(singular_values3, 'filled');
title('Singular Values of X with Doubled Number of Antennas');
xlabel('Index');
ylabel('Singular Value');
% saveas(gcf,'singulardoubleM.png')

% Scenario 3: Small angles between sources
theta2 = [-5, 5]'; % Smaller angles between sources
[X4, ~, ~] = gendata(M, N, Delta, theta2, f, SNR);
singular_values4 = svd(X4);
figure;
stem(singular_values4, 'filled');
title('Singular Values of X with Small Angles Between Sources');
xlabel('Index');
ylabel('Singular Value');
% saveas(gcf,'singularsmallangles.png')

% Scenario 4: Small frequency difference
f2 = [0.12, 0.15]'; % Smaller frequency difference
[X5, ~, ~] = gendata(M, N, Delta, theta, f2, SNR);
singular_values5 = svd(X5);
figure;
stem(singular_values5, 'filled');
title('Singular Values of X with Small Frequency Difference');
xlabel('Index');
ylabel('Singular Value');
% saveas(gcf,'singularsmallfreqs.png')
%% Estimation of directions with ESPRIT
estimated_theta = esprit(X, d);

%% Estimation of frequencies with ESPRIT
estimated_f = espritfreq(X, d);

%% Joint diagonalization
% Generate data with Khatri rao structure
% [X, ~, ~] = gendata(M, N, Delta, theta, f, SNR);
m = 3;
[estimated_theta, estimated_f] = joint(X,d,m);

%% Comparison ESPRIT and Joint Diagonalization
close all;
clc;
clear;

% Parameters
% Define parameters
M = 10; % Number of antennas
N = 20; % Number of samples
d = 2; % Number of sources
Delta = 0.5; % Antenna spacing in wavelengths
theta = [10, -15]'; % Directions of the sources in degrees
f = [0.3, 0.1]'; % Normalized frequencies of the sources
SNR_values = 0:4:20; % SNR values to test
num_runs = 1000; % Number of test runs

% Preallocate arrays to store results
theta_esprit = zeros(num_runs, length(SNR_values), 2);
freq_esprit = zeros(num_runs, length(SNR_values), 2);
theta_joint = zeros(num_runs, length(SNR_values), 2);
freq_joint = zeros(num_runs, length(SNR_values), 2);

% Loop over SNR values
for idx = 1:length(SNR_values)
    SNR = SNR_values(idx);
    for run = 1:num_runs
        % Generate data
        [X, ~, ~] = gendata(M, N, Delta, theta, f, SNR);

        % ESPRIT
        theta_esprit(run, idx, :) = esprit(X, 2);

        % espritfreq
        freq_esprit(run, idx, :) = espritfreq(X, 2);

        % Joint Estimation
        [theta_joint(run, idx, :), freq_joint(run, idx, :)] = joint(X, 2, 3);
    end
end

% Calculate mean and standard deviation for each algorithm
mean_theta_esprit = mean(theta_esprit, 1);
std_theta_esprit = std(theta_esprit, 0, 1);
mean_freq_esprit = mean(freq_esprit, 1);
std_freq_esprit = std(freq_esprit, 0, 1);
mean_theta_joint = mean(theta_joint, 1);
std_theta_joint = std(theta_joint, 0, 1);
mean_freq_joint = mean(freq_joint, 1);
std_freq_joint = std(freq_joint, 0, 1);

% Plotting
figure;
subplot(2, 1, 1);
errorbar(SNR_values, squeeze(mean_theta_esprit(:, :, 1)), squeeze(std_theta_esprit(:, :, 1)));
hold on;
errorbar(SNR_values, squeeze(mean_theta_joint(:, :, 1)), squeeze(std_theta_joint(:, :, 1)));
xlabel('SNR (dB)');
ylabel('Mean and STD of Estimated Angles (degrees)');
legend('ESPRIT', 'Joint');
title('Angle Estimation Performance');

subplot(2, 1, 2);
errorbar(SNR_values, squeeze(mean_freq_esprit(:, :, 1)), squeeze(std_freq_esprit(:, :, 1)));
hold on;
errorbar(SNR_values, squeeze(mean_freq_joint(:, :, 1)), squeeze(std_freq_joint(:, :, 1)));
xlabel('SNR (dB)');
ylabel('Mean and STD of Estimated Frequencies');
legend('espritfreq', 'Joint');
title('Frequency Estimation Performance');
%% Channel equalization
close all;
clc;
clear;

% Parameters
N = 500; % N QPSK symbolsN QPSK symbols
P = 4; % Sampling rate of P
sigma = .5; % std of complex Gaussian noise
m = 2; % stacking m times

% Generate random bits for in-phase and quadrature components
bits_I = randi([0, 1], 1, N);
bits_Q = randi([0, 1], 1, N);

% Map bits to QPSK symbols
s = (2 * bits_I - 1) / sqrt(2) + 1j * (2 * bits_Q - 1) / sqrt(2);
x = gendata_conv(s,P,N,sigma);

%Data matrix
X = zeros(m*P, N-1);
for i = 1:N-1
    X(:,i) = x(1+(i-1)*P:(i+1)*P);
end

% create H matrix
H = zeros(m*P,1+m-1);
if P == 4
    H(1:P, 2) = [1;-1;1;-1];
    H(P+1:2*P, 1) = [1;-1;1;-1];
else
    H(1:P, 2) = [1; 1;-1; -1;1; 1;-1;-1];
    H(P+1:2*P, 1) = [1; 1;-1; -1;1; 1;-1;-1];
end

% Zero Forcing
H_inv = pinv(H);
w_conj = H_inv(1,:);
est_s_zf = w_conj*X;

% Plot scatterplot in complex plane
real_parts = real(est_s_zf);
imaginary_parts = imag(est_s_zf);
figure;
scatter(real_parts, imaginary_parts, 'filled');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Scatter Plot of the Receivers');
grid on;
hold on;

% Wiener receiver
W = inv( H*H'+ sigma^2 * eye(size(H,1)) ) * H;
est_s_wf = W' * X;
est_s_wf = est_s_wf(1,:);

% Plot scatterplot in complex plane
real_parts = real(est_s_wf);
imaginary_parts = imag(est_s_wf);
scatter(real_parts, imaginary_parts, 'filled');

% Plot QPSK symbols
theta = [1/4, 3/4, 5/4, 7/4] * pi;
unit_circle = exp(1i * theta);
real_parts = real(unit_circle);
imaginary_parts = imag(unit_circle);
scatter(real_parts, imaginary_parts, 'filled', 'yellow');
legend({'ZF','WF', 'QPSK True points'});