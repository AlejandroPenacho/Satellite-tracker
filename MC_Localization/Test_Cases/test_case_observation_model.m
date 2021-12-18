function valid = test_case_observation_model()

    % load test data
    load test_case_obs.mat;

    % import global variables
    global map
    global M
    
    % set global variables
    map = W;
    M = size(S_BAR, 2);

    % error containers
    errs = ones(1,NUM_TEST);
    z= zeros(2,PARTICLE_NUM,NUM_TEST);
    
    % run test cases
    for i = 1 : NUM_TEST
        
        % get test data for current iteration
        s = S_BAR(:,:,i);
        
        % get measurement estimates
        z(:,:,i) = observation_model(s, J(i));
        
        % keep error statistics
        err = z(:,:,i) - Z(:,:,i);
        errs(i) = norm(err);
    end   
    
    % compute MSE
    mse_err = mean(errs.^2);
    THRESH_VALID = 1e-20;
    valid = mse_err < THRESH_VALID;

end