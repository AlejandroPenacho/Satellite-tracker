classdef ParticleFilter
    % Represents a particle filter, following a concrete satellite in Earth
    % orbit. 
    
    properties
        target              % The solution of an ODE, providing state of the target at a given time
        S                   % A DxN array with the state of the N particles, where the state has D dimensions
        n_particles         % The number of particles the filter is using
        time                % The current clock time of the filter, in seconds
        filter_state        % An struct with some variables that the filter uses to tune itself
        three_dimensional   % A boolean indicating whether 3- dimensions are bing used
        observation_system  % The observation system used by the filter, that is, how measurements are obtained
        predictor           % Updates the state of the particles using its own orbital motion model adding some uncertainty
        updater             % Uses the real measurements to resample the particle set, so its distributions approaches the real one
    end
    
    methods
        function obj = ParticleFilter(n_particles, IC, three_dimensional, dispersion_models, ground_stations, target)
            % Function used to create the filter

            obj.n_particles = n_particles;
            obj.three_dimensional = three_dimensional;

            % There are different methods to initialize the particles 
            obj.S = handle_IC(IC,n_particles);

            obj.time = 0;

            obj.target = target;

            % The predictor, observations system and the updater are
            % initialized

            obj.predictor = Predictor(n_particles, dispersion_models, three_dimensional);
            obj.observation_system = ObservationSystem(three_dimensional, ground_stations);
            obj.updater = Updater(n_particles, three_dimensional);

            % Initialization of the filter_state
            obj.filter_state = struct( ...
                "detected", false, ...                               % Whether the target was seen in last measurement
                "time_since_detection", 0, ...                      % Time since the last step with no detection
                "active_gs", true(length(ground_stations),1), ...   % Which ground stations are currently seeing the target
                "weight_variance", ones(1,obj.n_particles), ...     % Variance of the weights of the particles in the resampling
                "detection_status", "first_contact" ...
                );
        end

        function [state, covariance] = get_estimation(obj)

            % Provides the estimated state and the covariance of this
            % estimation of the target

            state = sum(obj.S,2)/obj.n_particles;
            delta = obj.S - repmat(state,1,obj.n_particles);
            covariance = delta*delta'/obj.n_particles;
        end
        
        function obj = step(obj, delta_t, skip_update, real_obs, real_detection)
            % One step of the filter, using prediction and update to obtain
            % a new set of particles tracking the target.
            
            % If the measurements of the target are not provided, they can
            % be obtained from the internal "target" variable, and using
            % the observational system of the particle.
            if nargin == 2
                obj.time = obj.time + delta_t;
                X = deval(obj.target, obj.time);

                [real_obs, real_detection] = obj.observation_system.get_measurements(X, obj.time);
                skip_update = false;
            end

            if skip_update
                obj.filter_state.time_since_detection = 0;
                obj.filter_state.detection_status = obj.update_detection_status();
                [obj.S, obj.filter_state] = obj.predictor.predict(obj.S, delta_t, obj.filter_state);
                return
            end

            % The ground stations actively detecting the target are
            % inserted in the filter state, in order to be used for other
            % operations

            obj.filter_state.active_gs = reshape(real_detection,1,[]) == 1;

            % Update the time since detection of the filter

            if sum(obj.filter_state.active_gs) > 0
                obj.filter_state.time_since_detection = obj.filter_state.time_since_detection + delta_t;
            else 
                obj.filter_state.time_since_detection = 0;
            end

            obj.filter_state.detection_status = obj.update_detection_status();

            real_obs = real_obs(:,:,real_detection==1);

            [S_bar, obj.filter_state] = obj.predictor.predict(obj.S, delta_t, obj.filter_state);


            [particle_obs, ~] = obj.observation_system.get_measurements(S_bar, obj.time, obj.filter_state.active_gs);
           
            R = obj.observation_system.get_R(obj.filter_state.active_gs);

            [obj.S, obj.filter_state] = obj.updater.update(S_bar, real_obs, particle_obs, R, obj.filter_state);
        end

        function detection_status = update_detection_status(obj)
            if obj.filter_state.time_since_detection > 20
                detection_status = "normal_contact";
            elseif obj.filter_state.time_since_detection > 0
                if obj.filter_state.detection_status == "first_contact"
                    detection_status = "first_contact";
                else
                    detection_status = "recovery";
                end
            else
                detection_status = "no_sight";
            end
                        
        end

        function plot_state(obj)
            if obj.three_dimensional
                % Earth
                
                [X,Y,Z] = sphere(20);
                surf(6371*X,6371*Y,6371*Z, "FaceColor",[0.6863, 0.3690, 0], ...
                    "FaceAlpha",0.5);

                hold on


                % Real target
                X = deval(obj.target, obj.time);
                scatter3(X(1), X(2),X(3), 20, [0.4940, 0.1840, 0.5560])

                % Stations and vision lines

                i_earth = 0;

                for i=1:obj.observation_system.n_gs
                    gs_lat = obj.observation_system.gs_location(1,i);
                    gs_long = obj.observation_system.gs_location(2,i) + 2*pi/(24*3600) * obj.time;

                    radial_vector = [
                        cos(gs_lat).*cos(gs_long);
                        cos(i_earth).*cos(gs_lat).*sin(gs_long) - sin(i_earth).*sin(gs_lat);
                        sin(i_earth).*cos(gs_lat).*sin(gs_long) + cos(i_earth).*sin(gs_lat)
                    ];


                    if obj.filter_state.active_gs(i)
                        color = [0.4660, 0.6740, 0.1880];
                    else
                        color = [0.6350, 0.0780, 0.1840];
                    end

                    scatter3(6371*radial_vector(1), 6371*radial_vector(2), 6371*radial_vector(3), 30, color, "filled");
                    plot3([6371*radial_vector(1), X(1)], ...
                          [6371*radial_vector(2), X(2)], ...
                          [6371*radial_vector(3), X(3)], ...
                          "color", color, "LineWidth",1.5);
                end

                % Particles
                scatter3(obj.S(1,:), obj.S(2,:), obj.S(3,:), 6, [0, 0.4470, 0.7410], "filled")

                daspect([1 1 1]);

                hold off                

            else
                % Earth
                earth = zeros(2,1001);
                for i=1:1001
                    earth(:,i) = 6371*[cos(i*2*pi/1000); sin(2*pi*i/1000)];
                end

                plot(earth(1,:), earth(2,:), "color", [0.6863, 0.3690, 0], "LineWidth", 1.5)

                hold on


                % Real target
                X = deval(obj.target, obj.time);
                scatter(X(1), X(2), 20, [0.4940, 0.1840, 0.5560])

                % Stations and vision lines
                for i=1:obj.observation_system.n_gs
                    gs_longitude = obj.observation_system.gs_location(i) + 2*pi/(24*3600) * obj.time;

                    if obj.filter_state.active_gs(i)
                        color = [0.4660, 0.6740, 0.1880];
                    else
                        color = [0.6350, 0.0780, 0.1840];
                    end

                    scatter(6371*cos(gs_longitude), 6371*sin(gs_longitude), 30, color, "filled");
                    plot([6371*cos(gs_longitude), X(1)], [6371*sin(gs_longitude), X(2)], "color", color);
                end

                % Particles
                scatter(obj.S(1,:), obj.S(2,:), 6, [0, 0.4470, 0.7410], "filled")

                daspect([1 1 1]);

                hold off
            end
        end
    end
end

