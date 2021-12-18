% This function is the implementation of the measurement model.
% The bearing should be in the interval [-pi,pi)
% Inputs:
%           S(t)                           4XM
%           j                              1X1
% Outputs:  
%           z_j                              2XM
function z_j = observation_model(S, j)

    global map % map including the coordinates of all landmarks | shape 2Xn for n landmarks
    global M % number of particles

    % YOUR IMPLEMENTATION

    z_j = [
        sqrt((S(1,:)- map(1,j)).^2 + (S(2,:) - map(2,j)).^2 );
        mod(atan2(map(2,j) - S(2,:), map(1,j) - S(1,:)) - S(3,:) + pi, 2 * pi) - pi
    ];



end
