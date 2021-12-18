function valid = test_case_weight()

    % load test data
    load test_case_weight.mat;
    
    % error containers
    errs = ones(1,NUM_TEST);
    S_bar_N = zeros(size(S_BAR));
    
    % run test cases
    for i = 1:NUM_TEST
        
        % get test data for current iteration
        s_bar = S_BAR(:,:,i);
        psi = reshape(PSI(:,:,i),1,NUM_OBS,PARTICLE_NUM);
        outlier = OUTLIERS(:,i)';
        
        % compute particle weights
        S_bar_N(:,:,i) = weight(s_bar,psi,outlier);
        
        % keep errror statistics
        errs(i) = norm(S_BAR_N(:,:,i) - S_bar_N(:,:,i));
    end
    
    % compute MSE
    mse_err = mean(errs.^2);
    THRESH_VALID = 1e-20;
    valid = mse_err < THRESH_VALID;
end