function x = gendata_conv(s,P,N,sigma)
    % s: signal of QPSK symbols
    % P: Sample at a rate 1/P
    % N: signal length
    % sigma: std of zero-mean complex Gaussian noise

    % Channel impulse response h(t) sampled at rate 1/P
    t = (0:(1/P):(1 - 1/P))';
    h = ones(size(t));
    h(t >= 0.25 & t < 0.5) = -1;
    h(t >= 0.5 & t < 0.75) = 1;
    h(t >= 0.75) = -1;

    % Convolution of the data symbols with the channel impulse response
    x = zeros(N*P,1);
    for i = 1:N
        x(1+(i-1)*P:i*P) = h * s(i);
    end

    % Additive Gaussian noise
    noise = sigma * (randn(size(x)) + 1i * randn(size(x))) / sqrt(2);

    % Received signal
    x = x + noise;
end