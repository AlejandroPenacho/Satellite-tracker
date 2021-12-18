function test_cases_lab2()
clc;
done = 1;

% all relevant functions
names = {'calculate odometry', 'predict','observation model','associate','weight'};
functions = {'calculate_odometry', 'predict','observation_model','associate','weight'};

% run all test cases
for i = 1 : length(functions)
    if eval(sprintf('test_case_%s()',functions{i}))
        fprintf('Your %s function seems to be fine! \n',names{i});
    else
        fprintf('Your %s function seems to be wrong, not running other test cases! \n',names{i});
        done = 0;
        break;
    end
end
if done
    disp('Congratulations! It seems that you are on a good track!');
else
    disp('Oops! It seems you still need some debugging!');
end
end