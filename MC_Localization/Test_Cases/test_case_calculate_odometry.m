function valid =test_case_calculate_odometry()

    % load test data
    load test_case_odo.mat;
    
    % error containers
    errs = ones(1,NUM_TEST);
    v = zeros(1,NUM_TEST);
    omega= zeros(1,NUM_TEST);
    
    % run test cases
    for i = 1 : NUM_TEST
        
        % compute odometry
        [v(i),omega(i)] = calculate_odometry(E_R(i),E_L(i),E_T,B,R_R,R_L,delta_t);
        
        % keep error statistics
        errs(i) = norm([v(i);omega(i)] - [V(i);OMEGA(i)]);
    end
    
    % compute MSE
    mse_err = mean(errs .^2);
    THRESH_VALID = 1e-20;
    valid = mse_err < THRESH_VALID;
end