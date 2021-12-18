function valid =test_case_predict()

    % load test data
    load test_case_prediction.mat;

    % import global variables
    global M
    global R

    % set global variables
    M = PARTICLE_NUM;

    % error containers
    errs = ones(1,NUM_TEST);
    S_bar = zeros(4,PARTICLE_NUM,NUM_TEST);

    % run test cases
    for i = 1 : NUM_TEST

        % get test data for current iteration
        s = S(:,:,i);
        v = V(i);
        R = diag(RT(:,i));
        omega = OMEGA(i);

        % get predictions
        for t = 1 : RUN_TESTS_AVERAGE
            [s_bar] = predict(s, v, omega, delta_t);
            S_bar(:,:,i) = S_bar(:,:,i) + s_bar;
        end

        % keep error statistics
        S_bar(:,:,i) = S_bar(:,:,i) ./ RUN_TESTS_AVERAGE;
        err = S_BAR(1:3,:,i) -S_bar(1:3,:,i);
        errs(i) = mean(sqrt(sum((err).^2,1)));
    end

    % compute MSE
    mse_err = mean(errs .^2);
    THRESH_VALID = 1e-3;
    valid = mse_err < THRESH_VALID;
end