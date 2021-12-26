classdef ObservationSystem
    %OBSERVATIONSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ground_stations
        n_gs
        three_dimensional
        std_R
    end
    
    methods
        function obj = ObservationSystem(three_dimensional, ground_stations)
            % ground_station: a cell for each station, defined by position
            % [latitude; longitude] and precision [distance; alpha; beta].
            % For 2D, there is no latitude nor alpha angle

            obj.n_gs = length(ground_stations);
            obj.ground_stations = ground_stations;
            obj.three_dimensional = three_dimensional;

            obj.std_R = [];

            for i=1:n_gs
                obj.std_R = blkdiag(obj.std_R, 1/12 * ground_stations{i}.precision.^2);
                obj.ground_stations{i}.rounding = floor(-log10(ground_stations{i}.precision));
            end
        end
        
        function [measurements, detection] = get_observations(obj, particles, time)
            
            if obj.three_dimensional
                [measurements, detection] = get_3D_observations(obj, particles, time);
            else
                [measurements, detection] = get_2D_observations(obj, particles, time);
            end
        end

        function [measurements, detection] = get_2D_observations(obj, particles, time)

            n_particles = size(particles, 2);
        
            earth_radius = 6371;
        
            gs_long = obj.ground_stations(1,:) + 2*pi/(24*3600) * time;
        
            radial_vector = [
                cos(gs_long);
                sin(gs_long)
            ];
        
            az_vector = [
                -sin(gs_long);
                cos(gs_long)
            ];
        
            gs_position = earth_radius * radial_vector;
        
            delta_x = repmat(particles(1:2,:),1,1, obj.n_gs) - ...
                      repmat(reshape(gs_position,2,1, obj.n_gs),1,n_particles,1);
        
        
            distance = vecnorm(delta_x);
        
            radial_components = pagemtimes(reshape(radial_vector,1,2, obj.n_gs), delta_x./distance);
            azimuthal_components = pagemtimes(reshape(az_vector,1,2, obj.n_gs), delta_x./distance);
        
            beta = asin(azimuthal_components);
        
            measurements = [distance; beta];
        
            detection = radial_components >= 0;

        end
    end
end

