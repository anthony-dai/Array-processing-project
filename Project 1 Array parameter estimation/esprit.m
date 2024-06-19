function theta = esprit(X, d)
%     % X - received signal matrix (MxN)
%     % d - number of sources

    % Partition the signal subspace into two subarrays
    X_s = X(1:end-1, :);
    Y_s = X(2:end, :);
    Z = [X_s ; Y_s];
    
    % Compute SVD
    [U,S,~] = svd(Z);
    [~, idx] = sort(diag(S), 'descend');
    U = U(:, idx);
    U_hat = U(:,1:d);
%     disp(size(X))
%     disp(size(U_hat))
    
    % Split U hat in half
    middle = size(U_hat,1)/2;
    Ux_hat = U_hat(1:middle, :);
    Uy_hat = U_hat((middle+1):end, :);
    
    Ux_crossUy = pinv(Ux_hat) * Uy_hat;
    
    % eigenvalue decomposition
    [~, D] = eig(Ux_crossUy);
    %phis = sort(diag(D), 'descend');
    phis = diag(D);
    tmp1 = log(phis);
    tmp2 = -1j*tmp1/(2 * pi * 0.5);
    theta = real(asin(tmp2)) * 180 / pi;
end