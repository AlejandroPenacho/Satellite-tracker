classdef MultitargetFilter
    %MULTITARGETSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time
        targets
        three_dimensional
        std_filter
        ground_stations
        observation_system
        current_filters
    end
    
    methods
        function obj = MultitargetFilter(n_particles, three_dimensional, dispersion_models, ground_stations, targets)
            obj.std_filter = struct( ...
                "n_particles", n_particles, ...
                "three_dimensional", three_dimensional, ...
                "dispersion_models", dispersion_models);
            
            obj.ground_stations = ground_stations;
            obj.targets = targets;
            obj.three_dimensional = three_dimensional;
            obj.observation_system = ObservationSystem(three_dimensional, ground_stations);
            obj.time = 0;

        end
        
        function obj = step(obj, delta_t)

            p_threshold = 0.00001;

            obj.time = obj.time + delta_t;

            [observations, detection_array] = obj.generate_real_observations();
            [target_measurements, sigmas] = obj.generate_target_measurements();

            filter_ref = 1:length(obj.current_filters);

            for i = 1:size(observations,2)
                if sum(sum(detection_array(:,i,:))) == 0
                    continue
                end

                [best_p_index, best_p] = obj.get_best_fit(observations(:,i,:), detection_array(:,i,:), target_measurements, sigmas, filter_ref);
                if best_p > p_threshold && ~isempty(filter_ref)
                    filter_index = filter_ref(best_p_index);
                    obj.current_filters{filter_index} = obj.current_filters{filter_index}.step(delta_t, false, observations(:,i,:), detection_array(:,i,:));
                    filter_ref(best_p_index) = [];
                else
                    obj.current_filters{end+1} = obj.create_new_filter(observations(:,i,:), detection_array(:,i,:));
                end
            end

            for j=1:length(filter_ref)
                obj.current_filters{filter_ref(j)} = obj.current_filters{filter_ref(j)}.step(delta_t, true);
            end

        end

        function filter = create_new_filter(obj, observation, detection_array)

            detected = false;
            for gs_index=1:size(detection_array,3)
                if detection_array(1,1,gs_index)
                    detected = true;
                    break
                end
            end
            
            if ~detected
                error("Damn it")
            end

            observation = observation(:,:,gs_index);
            d = observation(1);
            v = sqrt(398600/d);
            
            if obj.three_dimensional
                [gs_lat, gs_long] = obj.ground_stations.location(:, gs_index);
                gs_long = gs_long + obj.time*(2*pi/24*3600);
                i_earth = 0.4;
                T = obj.observation_system.obtain_transform_matrix(gs_long, gs_lat, i_earth);

                alpha = observation(2);
                beta = observation(3);

                delta_x_radial_coord = d*[cos(alpha)*cos(beta); cos(alpha)*sin(beta); sin(alpha)];
                delta_x = (T')*delta_x_radial_coord;
                x_gs = d*T(:,1);
                x_sat = x_gs + delta_x;

                IC = struct("X", x_sat, "mode", "3d-start");
                filter = ParticleFilter(obj.n_particles, IC, true, obj.dispersion_models, obj.ground_stations, false);
            
            else
                gs_long = obj.ground_stations{gs_index}.location;
                gs_long = gs_long + obj.time*(2*pi/(24*3600));
                T = obj.observation_system.obtain_transform_matrix(gs_long);
                beta = observation(2);

                delta_x_radial_coord = d*[cos(beta); sin(beta)];
                delta_x = (T)*delta_x_radial_coord;
                x_gs = 6371*T(:,1);
                x_sat = x_gs + delta_x;

                IC = struct("X", x_sat, ...
                            "mode", "2d-start");
                filter = ParticleFilter(obj.std_filter.n_particles, IC, false, obj.std_filter.dispersion_models, obj.ground_stations, false);                
            end

            filter.time = obj.time;
        end

        function [observations, detection_array] = generate_real_observations(obj)

            particles = [];

            for i=1:length(obj.targets)
                if obj.time >= obj.targets{i}.detection_time
                    particles = [particles, deval(obj.targets{i}.solved_ode, obj.time)];
                end
            end

            [observations, detection_array] = obj.observation_system.get_measurements(particles, obj.time);

        end

        function [measurements, sigmas] = generate_target_measurements(obj)

            particles = [];
            sigmas = cell(length(obj.current_filters),1);

            for i=1:length(obj.current_filters)
                [mu, sigma] = obj.current_filters{i}.get_estimation();
                particles = [particles, mu];

                C = obj.observation_system.get_C(mu);

                sigmas{i} = sigma;
            end

            if ~isempty(obj.current_filters)
                [measurements, ~] = obj.observation_system.get_measurements(particles, obj.time);
            else
                measurements = [];
            end
        end        

        function [best_index, best_psi] = get_best_fit(obj, observation, detection_array, target_measurements, sigmas, filter_ref)


            if isempty(target_measurements)
                best_index = 0;
                best_psi = 0;
                return
            end
            
            % best_index = 1;
            % best_psi = 1;
            % return

            n_targets = length(sigmas);

            target_measurements = target_measurements(:,:,detection_array);
            observation = observation(:,:,detection_array);
    
            z = reshape(permute(target_measurements,[1,3,2]),[],n_targets) - repmat(reshape(observation,[],1),1,n_targets);

            best_psi = 0;
            best_index = 0;

            f_sigma = diag([50, 0.1]);
            % f_sigma = diag([0.01, 0.01]);

            for i=1:length(filter_ref)
                filter_index = filter_ref(i);
                % psi = 1/sqrt(2*pi*det(sigmas{filter_index}))* ...
                %     exp(-0.5*z(:,filter_index)'*(sigmas{filter_index}^-1)*z(:,filter_index));
                                            % Substitute sigma by
                                            % C*sigma*C'

                psi = 1/sqrt(2*pi*det(f_sigma))* ...
                    exp(-0.5*z(:,filter_index)'*(f_sigma^-1)*z(:,filter_index));

                if psi > best_psi
                    best_psi = psi;
                    best_index = i;
                end
            end

        end

        function plot(obj)
            if obj.three_dimensional
                error("No plot in 3-D")
            end

            % Earth

            earth = zeros(2,1001);
            for i=1:1001
                earth(:,i) = 6371*[cos(i*2*pi/1000); sin(2*pi*i/1000)];
            end

            plot(earth(1,:), earth(2,:), "color", [0.6863, 0.3690, 0], "LineWidth", 1.5)

            hold on

            % Real objects
            for i=1:length(obj.targets)
                if obj.time >= obj.targets{i}.detection_time
                    X = deval(obj.targets{i}.solved_ode, obj.time);
                    scatter(X(1), X(2), 20, [0.4940, 0.1840, 0.5560])
                end
            end

            % Estimations

            for i=1:length(obj.current_filters)
                scatter(obj.current_filters{i}.S(1,:), obj.current_filters{i}.S(2,:), 6, [0, 0.4470, 0.7410], "filled")
            end

            hold off

        end

    end
end

