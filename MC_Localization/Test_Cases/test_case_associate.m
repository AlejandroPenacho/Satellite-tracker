function valid = test_case_associate()

    % load test data
    load test_case_associate.mat;

    % import global variables
    global lambda_psi % threshold on average likelihood for outlier detection
    global Q % covariance matrix of the measurement model
    global M % number of particles
    global N % number of landmarks
    global map;
    map=W;
    % set global variables
    M = size(S_BAR, 2);
    N = size(W, 2);

    % error containers 
    errs_psi = ones(1, NUM_TEST);
    errs_outliers = ones(1, NUM_TEST);
    psi = zeros(NUM_OBS, PARTICLE_NUM,NUM_TEST);
    outliers = zeros(NUM_OBS, NUM_TEST);
    
    % run test cases
    for i = 1 : NUM_TEST
        
        % get test data for current iteration
        s_bar = S_BAR(:,:,i);
        Q = diag(QT(:,i));
        z_i = Z(:,:,i);
        lambda_psi = LAMBDA(i);
        
        % perform association
        [outliers(:,i),psi(:,:,i), ~] = associate(s_bar,z_i);
        
        % keep error statistics
        err_outlier = OUTLIERS(:,i) - outliers(:,i);
        err_psi = PSI(:,:,i) - psi(:,:,i);
        errs_psi(i) = norm(err_psi);
        errs_outliers(i) = norm(err_outlier);
    end   
    
% check MSE
mse_err_outliers = mean(errs_outliers.^2);
mse_err_psi = mean(errs_psi.^2);
THRESH_VALID = 1e-20;
valid_outliers = mse_err_outliers < THRESH_VALID;
if ~ valid_outliers
    fprintf('outliers computed in associate.m seem to be wrong, mse=%f\n', mse_err_outliers);
end
valid_psi = mse_err_psi < THRESH_VALID;
if ~ valid_psi
    fprintf('psi s computed in associate.m seem to be wrong, mse=%f\n',mse_err_psi);
end

% associate function valid if psi and outliers correct
valid = valid_outliers && valid_psi;
end
