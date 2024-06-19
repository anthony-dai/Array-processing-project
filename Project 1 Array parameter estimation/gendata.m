function [X,A,S] = gendata(M,N,Delta,theta,f,SNR)
    % M - number of antennas
    % N - number of samples
    % Delta - antenna spacing (in wavelengths)
    % theta - directions of the sources (in degrees)
    % f - normalized frequencies of the sources
    % SNR - signal-to-noise ratio per source

    d = length(theta);  % Number of sources
    theta = theta(:)';  % Ensure theta is a row vector
    f = f(:)';          % Ensure f is a row vector

    % Convert angles from degrees to radians
    theta_rad = theta * pi / 180;
    
    % Generate steering matrix A
    A = zeros(M, d);
    for i = 1:d
        A(:, i) = exp(1j * 2 * pi * Delta * (0:M-1)' * sin(theta_rad(i)));
    end
    
    % Generate source signals S
    S = zeros(d, N);
    for i = 1:d
        S(i, :) = exp(1j * 2 * pi * f(i) * (0:N-1));
    end
    
    % Compute the noise power from SNR
    signal_power = 1; % Assume signal power is 1
    noise_power = signal_power / (10^(SNR / 10));
    
    % Generate noise
    noise = sqrt(noise_power/2) * (randn(M, N) + 1j * randn(M, N));

    % Generate received signal matrix X
    X = A * S + noise;
end