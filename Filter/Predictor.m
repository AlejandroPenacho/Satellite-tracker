classdef Predictor
    % This class corresponds to the prediction phase of the filter.
    % It requires only one method: predict, that takes the particles and
    % updates their positions according to its model.
    
    %TODO: Q should scale with timestep!!!!

    properties
        Q
        n_particles
        three_dimensional
        mu_earth
    end
    
    methods
        function obj = Predictor(n_particles, Q, three_dimensional)
            %PREDICTOR Construct an instance of this class
            %   Detailed explanation goes here

            obj.mu_earth = 398600;
            obj.three_dimensional = three_dimensional;
            obj.n_particles = n_particles;

            if three_dimensional
                assert(size(Q,1) == 6 && size(Q,2) == 6)
            else
                assert(size(Q,1) == 4 && size(Q,2) == 4)
            end

            obj.Q = Q;
        end
        
        function S_bar = predict(obj, S, delta_t)
            % Given the current set of particles and a delta_t, update the 
            % state of the particles

            if ~obj.three_dimensional
                force_factor = obj.mu_earth ./ (S(1,:).^2+S(2,:).^2).^(3/2);
    
                S_dot = [
                    S(3,:);
                    S(4,:);
                    - force_factor .* S(1,:);
                    - force_factor .* S(2,:)
                ];

                S_bar = S + S_dot*delta_t + [
                    mvnrnd([0;0;0;0], obj.Q, obj.n_particles)'
                ];

            else

                force_factor = obj.mu_earth ./ (S(1,:).^2 + S(2,:).^2 + S(3,:).^2).^(3/2);
                
                S_dot = [
                    S(4,:);
                    S(5,:);
                    S(6,:);
                    - force_factor .* S(1,:);
                    - force_factor .* S(2,:);
                    - force_factor .* S(3,:)
                ];

                S_bar = S + S_dot*delta_t + [
                    mvnrnd([0;0;0;0;0;0], obj.Q, obj.n_particles)'
                ];

            end

        end
    end
end

