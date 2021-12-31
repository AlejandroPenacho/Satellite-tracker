classdef Predictor
    % This class corresponds to the prediction phase of the filter.
    % It requires only one method: predict, that takes the particles and
    % updates their positions according to its model.
    
    %TODO: Q should scale with timestep!!!!

    properties
        dispersion_models
        n_particles
        three_dimensional
        mu_earth
    end
    
    methods
        function obj = Predictor(n_particles, dispersion_models, three_dimensional)
            %PREDICTOR Construct an instance of this class
            %   Detailed explanation goes here

            obj.mu_earth = 398600;
            obj.three_dimensional = three_dimensional;
            obj.n_particles = n_particles;

            % if three_dimensional
            %     assert(size(Q,1) == 6 && size(Q,2) == 6)
            % else
            %     assert(size(Q,1) == 4 && size(Q,2) == 4)
            % end

            obj.dispersion_models = dispersion_models;
        end
        
        function [S_bar, filter_state] = predict(obj, S, delta_t, filter_state)
            % Given the current set of particles and a delta_t, update the 
            % state of the particles


            switch filter_state.detection_status

                case "first_contact"
                    Q = obj.dispersion_models.first_contact;

                case "no_sight"
                    Q = obj.dispersion_models.no_sight;

                case "recovery"
                    Q = obj.dispersion_models.recovery;

                case "normal_contact"
                    Q = obj.dispersion_models.standard;
                otherwise
                    error("Invalid detection status")
            end


            if ~obj.three_dimensional
                force_factor = obj.mu_earth ./ (S(1,:).^2+S(2,:).^2).^(3/2);
    
                S_dot = [
                    S(3,:);
                    S(4,:);
                    - force_factor .* S(1,:);
                    - force_factor .* S(2,:)
                ];

                S_bar = S + S_dot*delta_t + [
                    mvnrnd([0;0;0;0], Q, obj.n_particles)'
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
                    mvnrnd([0;0;0;0;0;0], Q, obj.n_particles)'
                ];

            end

        end
    end
end

