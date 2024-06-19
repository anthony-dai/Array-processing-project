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
M = 3; % Number of antennas
N = 20; % Number of samples
Delta = 0.5; % Antenna spacing in wavelengths
theta = [-20, 30]'; % Directions of the sources in degrees
f = [0.1, 0.12]'; % Normalized frequencies of the sources
SNR_values = 0:4:20; % SNR values to test
num_runs = 1000; % Number of test runs

% Preallocate arrays to store results
theta_esprit = zeros(num_runs, 2, length(SNR_values));
freq_esprit = zeros(num_runs, 2, length(SNR_values));
theta_joint = zeros(num_runs, 2, length(SNR_values));
freq_joint = zeros(num_runs, 2, length(SNR_values));

% Loop over SNR values
for idx = 1:length(SNR_values)
    SNR = SNR_values(idx);
    for run = 1:num_runs
        % Generate data
        [X, ~, ~] = gendata(M, N, Delta, theta, f, SNR);

        % ESPRIT for angles
        theta_est = esprit(X, 2);
        theta_esprit(run, :, idx) = theta_est(:);

        % ESPRIT for frequencies
        freq_est = espritfreq(X, 2);
        freq_esprit(run, :, idx) = freq_est(:);

        % Joint Estimation
        [theta_joint_est, freq_joint_est] = joint(X, 2, 3);
        theta_joint(run, :, idx) = theta_joint_est(:);
        freq_joint(run, :, idx) = freq_joint_est(:);
    end
end

% Average results over runs for each SNR value
avg_theta_esprit = squeeze(mean(theta_esprit, 1));
avg_freq_esprit = squeeze(mean(freq_esprit, 1));
avg_theta_joint = squeeze(mean(theta_joint, 1));
avg_freq_joint = squeeze(mean(freq_joint, 1));

% Plot results
figure;
subplot(2,1,1);
plot(SNR_values, avg_theta_esprit(1, :), 'r', SNR_values, avg_theta_esprit(2, :), 'r--');
hold on;
plot(SNR_values, avg_theta_joint(1, :), 'b', SNR_values, avg_theta_joint(2, :), 'b--');
xlabel('SNR (dB)');
ylabel('Angle (degrees)');
legend('ESPRIT Angle 1', 'ESPRIT Angle 2', 'Joint Angle 1', 'Joint Angle 2');
title('Angle Estimation vs SNR');

subplot(2,1,2);
plot(SNR_values, avg_freq_esprit(1, :), 'r', SNR_values, avg_freq_esprit(2, :), 'r--');
hold on;
plot(SNR_values, avg_freq_joint(1, :), 'b', SNR_values, avg_freq_joint(2, :), 'b--');
xlabel('SNR (dB)');
ylabel('Frequency');
legend('ESPRIT Freq 1', 'ESPRIT Freq 2', 'Joint Freq 1', 'Joint Freq 2');
title('Frequency Estimation vs SNR');


%%

% Generate data without noise
%[X, A, S] = gendata(M, N, Delta, theta, f, 1000); % High SNR for no noise
d = 2;
[X, A, S] = gendata(M, N, Delta, theta, f, 10); % Low SNR for no noise
%% Beamforming based on direction estimates
% Estimate angles using ESPRIT
theta_est_esprit = esprit(X, 2);

% Construct steering matrix
A_beam = zeros(M, d);
for i = 1:d
    phiZF_Angle = exp(1j*2*pi*Delta*sind(theta_est_esprit(i)));
    for k = 0:M-1
        A_beam(k+1, i) = phiZF_Angle^k;
    end
end

% Compute beamforming weights
W_esprit = (pinv(A_beam)).';

% Recover sources
S_esprit = W_esprit.' * X;

%% Beamforming based on frequency estimates
freq_est_freqesprit = espritfreq(X, 2);

% Construct steering matrix S_hat based on frequency estimates
S_hat = zeros(d, N); 
for i = 1:d
    phiZF_Freq = exp(1j*2*pi*freq_est_freqesprit(i));
    for k = 0:N-1
        S_hat(i, k+1) = phiZF_Freq^k;
    end
end

A_hat = X * pinv(S_hat);
W_espritfreq = pinv(A_hat);
S_espritfreq = W_espritfreq * X;
%%
% Compute spatial responses
theta_range = -90:1:90; % Angle range for plotting
spatial_response_freq = zeros(2, length(theta_range));
spatial_response_angle = zeros(2, length(theta_range));

for idx = 1:length(theta_range)
    theta = theta_range(idx);
    a_theta = exp(1j*2*pi*Delta*(0:M-1)'*sind(theta)); % Steering vector

    % Spatial response for angle-based beamformer
    spatial_response_angle(1,idx) = norm(W_esprit(:,1)' * a_theta,2);
    spatial_response_angle(2,idx) = norm(W_esprit(:,2)' * a_theta,2);
    % Spatial response for frequency-based beamformer
    spatial_response_freq(1,idx) = norm(W_espritfreq(2,:) * a_theta,2);
    spatial_response_freq(2,idx) = norm(W_espritfreq(1,:) * a_theta,2);
end

% Plot the spatial responses
figure;
plot(theta_range, 20*log10(spatial_response_freq(1,:)), 'r', 'LineWidth', 2);
hold on;
plot(theta_range, 20*log10(spatial_response_angle(1,:)), 'b', 'LineWidth', 2);
xlabel('Angle (degrees)');
ylabel('Spatial Response (dB)');
title('Spatial Responses of Zero-Forcing Beamformers');
legend('Frequency-based Beamformer', 'Angle-based Beamformer');
grid on;

% Plot the spatial responses
figure;
plot(theta_range, 20*log10(spatial_response_freq(2,:)), 'r', 'LineWidth', 2);
hold on;
plot(theta_range, 20*log10(spatial_response_angle(2,:)), 'b', 'LineWidth', 2);
xlabel('Angle (degrees)');
ylabel('Spatial Response (dB)');
title('Spatial Responses of Zero-Forcing Beamformers');
legend('Frequency-based Beamformer', 'Angle-based Beamformer');
grid on;


% Display the plots
hold off;
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