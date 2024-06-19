function f = espritfreq(X, d)
    % X - received signal matrix (MxN)
    % d - number of sources

    % Partition the signal subspace into two subarrays
    % Step 1: Form the data matrix Z
    X = X.';
%     Z = zeros(M, N - M);
%     for i = 1:M
%         Z(i, :) = X(i:N-M+i-1);
%     end

    X_s = X(1:end-1, :);
    Y_s = X(2:end, :);
    Z = [X_s ; Y_s];
    
    % Compute SVD
    [U,S,~] = svd(Z);
    [~, idx] = sort(diag(S), 'descend');
    U = U(:, idx);
    U_hat = U(:,1:d);
    
    % split in middle
    middle = size(U_hat,1)/2;
    Ux_hat = U_hat(1:middle, :);
    Uy_hat = U_hat((middle+1):end, :);
    Ux_crossUy = pinv(Ux_hat) * Uy_hat;
    
    % eigenvalue decomposition
    [~, D] = eig(Ux_crossUy); % V = T inv
    
%     [~, idx] = sort(diag(D), 'descend');
    phis = diag(D);
    omegas = angle(phis);
    for i = 1:d
        if omegas(i) < 0
            omegas(i) = omegas(i) + 2*pi;
        end
    end
    f = omegas/(2*pi);
end