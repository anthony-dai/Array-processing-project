function [theta,f] = joint(X,d,m)
    % X - received signal matrix (MxN)
    % d - number of sources
    % m - smoothing factor
    
    [M,N] = size(X);
    
    % Smoothening with factor m
    Z = zeros(M*m,N-m+1);
    for i=1:N-m+1
        for k=0:m-1
            Z(k*M+1:(k+1)*M,i)=X(:,i+k);
        end
    end
    
    % SVD on Z
    [U,S,~] = svd(Z);
    [~, idx] = sort(diag(S), 'descend');
    U = U(:, idx);
    U = U(:,1:d);
    
    % Splitting
    U_phix = U(1:M*(m-1),:);
    U_phiy = U(M+1:M*m,:);
    M_phi = pinv(U_phix)*U_phiy;
    U_thex = [];
    U_they = [];
    
    % 
    for i=0:m-1
        U_thex =  vertcat(U_thex,U(i*M+1:(i+1)*M-1,:));
        U_they =  vertcat(U_they,U(i*M+2:(i+1)*M,:));  
    end
    
    M_the = pinv(U_thex)*U_they;
    M = [M_the,M_phi];
    [~,D] = joint_diag(M,1.0e-8);
    [~,p] = size(D);
    Theta = D(:,1:p/2);
    Phi = D(:,p/2+1:p);
    
    % Return angles
    theta = asind(angle(eig(Theta))/pi);
    
    % Return frequencies
    omegas = angle(eig(Phi));
    for i = 1:d
        if omegas(i) < 0
            omegas(i) = omegas(i) + 2*pi;
        end
    end
    f = omegas/(2*pi);
end