% This function performs the prediction step.
% Inputs:
%           S(t-1)            4XN
%           v                 1X1
%           omega             1X1
% Outputs:   
%           S_bar(t)          4XN
function [S_bar] = predict(S, v, omega, delta_t)

    % Comment out any S_bar(3, :) = mod(S_bar(3,:)+pi,2*pi) - pi before
    % running the test script as it will cause a conflict with the test
    % function. If your function passes, uncomment again for the
    % simulation.

    global R % covariance matrix of motion model | shape 3X3
    global M % number of particles
    
    % YOUR IMPLEMENTATION

    
%    S_bar = S + [
%        v*delta_t*cos(S(3,:)) + normrnd(0, sqrt(R(1,1)), 1,M);
%        v*delta_t*sin(S(3,:)) + normrnd(0, sqrt(R(2,2)), 1,M);
%        omega*delta_t*ones(1,M) + normrnd(0, sqrt(R(3,3)), 1,M);
%        zeros(1,M)
%    ]; 

   S_bar = S + [
       v*delta_t*cos(S(3,:));
       v*delta_t*sin(S(3,:));
       omega*delta_t*ones(1,M);
       zeros(1,M)
   ] + [
       mvnrnd([0;0;0], R, M)';
       zeros(1,M)
   ];


    
    S_bar(3, :) = mod(S_bar(3,:)+pi,2*pi) - pi;
    
end