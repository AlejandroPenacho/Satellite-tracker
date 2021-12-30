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
        function obj = MultitargetSystem(n_particles, three_dimensional, dispersion_models, ground_stations, targets)
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
        
        function outputArg = step(obj, all_measurements)


            
            [filter_measurements, ~] = obj.observation_system.get_measurements

            for i=1:length(measurements)
                [measurement, detection] = obj.observation_system.get_measurement(all_measurements{i});

            end
        end

        function [best_fit, distance] = get_best_fit(obj,remaining_filters, measurement, detection)

            

            for i=1:length(remaining_filters)
                [state, ~] = remaining_filters.get_estimation();
                filter_obs = obj.remaining_filters.observation_system.get_measurements(state, obj.time, detection);
                reshaped_filter_obs = reshape(real_obs,[],1);
                reshaped_
            end
        end
    end
end

