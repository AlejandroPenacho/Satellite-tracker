% This function performs the ML data association
%           S_bar(t)                 4XM
%           z(t)                     2Xn
%           association_ground_truth 1Xn | ground truth landmark ID for
%           every measurement  1 (actually not used)
% Outputs: 
%           outlier                  1Xn    (1 means outlier, 0 means not outlier) 
%           Psi(t)                   1XnXM
%           c                        1xnxM (actually not ever used)
function [outlier, Psi, c] = associate(S_bar, z, association_ground_truth)
    if nargin < 3
        association_ground_truth = [];
    end

    global lambda_psi % threshold on average likelihood for outlier detection
    global Q % covariance matrix of the measurement model
    global M % number of particles
    global N % number of landmarks
   
    % YOUR IMPLEMENTATION
    
    Q_inv = Q^(-1);
    Q_inv_d = Q_inv(1);
    Q_inv_b = Q_inv(4);
    
    n_measurements = size(z,2);
    
    c = zeros(n_measurements,M);
        
    big =repmat(reshape(z',1,n_measurements,2,1),M,1,1,N);

    for i=1:N
        big(:,:,:,i) = big(:,:,:,i) - repmat(reshape(observation_model(S_bar,i)',M,1,2,1), 1, n_measurements, 1, 1);
    end
     
    big(:,:,2,:) =  mod(big(:,:,2,:) + pi, 2 * pi) - pi;
    
    big = (1/(2*pi*sqrt(det(Q))))*exp(-0.5*(Q_inv_d*big(:,:,1,:).^2+Q_inv_b*big(:,:,2,:).^2));
    
    [Psi, c] = max(big,[],4);
    
    
%     for p=1:M
%         for m=1:n_measures
%             psi = zeros(1,N);
%             
%             innovations = z(:,m) - observation_model(S_bar(:,p),lmk);
%             
%             for lmk=1:N
%                 innovation = z(:,m) - observation_model(S_bar(:,p),lmk);
%                 
%                 innovation(2) = mod(innovation(2) + pi, 2 * pi) - pi;
%                 
%                 psi(lmk) = (1/(2*pi*sqrt(det(Q))))*exp(-0.5*innovation'*(Q^-1)*innovation);
%             end
%             [Psi(m,p), c(m,p)] = max(psi);
%         end
%     end
    
    outlier = (sum(Psi,1)/M <lambda_psi);

    Psi = reshape(Psi',1,n_measurements,M);
    c = c';    
    
end
