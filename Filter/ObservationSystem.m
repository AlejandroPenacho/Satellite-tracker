classdef ObservationSystem
    %OBSERVATIONSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_gs
        three_dimensional
        std_R
        rounding
        gs_location
    end
    
    methods
        function obj = ObservationSystem(three_dimensional, ground_stations)
            % ground_station: a cell for each station, defined by position
            % [latitude; longitude] and precision [distance; alpha; beta].
            % For 2D, there is no latitude nor alpha angle

            obj.n_gs = length(ground_stations);
            obj.three_dimensional = three_dimensional;

            obj.std_R = {};
            obj.rounding = [];
            obj.gs_location = [];

            for i=1:obj.n_gs
                obj.gs_location = [obj.gs_location, ground_stations{i}.location];
                obj.rounding = [obj.rounding, floor(-log10(ground_stations{i}.precision))];
                obj.std_R{i} = diag(1/12 * ground_stations{i}.precision .^2);
            end
        end

        function R = get_R(obj, active_gs)
    
            R = [];
            for i=1:obj.n_gs
                if active_gs(i)
                    R = blkdiag(R, obj.std_R{i});
                end
            end
        end
        
        function [measurements, detection] = get_measurements(obj, particles, time, active_gs)
            
            if nargin == 3
                active_gs = true(obj.n_gs,1);
            end

            if obj.three_dimensional
                [measurements, detection] = obj.get_3D_observations(particles, time, active_gs);
            else
                [measurements, detection] = obj.get_2D_observations(particles, time, active_gs);
            end
        end


        function [measurements, detection] = get_2D_observations(obj, particles, time, active_gs)


            n_particles = size(particles, 2);

            n_active_gs = sum(active_gs);
        
            earth_radius = 6371;
        
            gs_long = obj.gs_location(active_gs) + 2*pi/(24*3600) * time;
        
            radial_vector = [
                cos(gs_long);
                sin(gs_long)
            ];
        
            az_vector = [
                -sin(gs_long);
                cos(gs_long)
            ];
        
            gs_position = earth_radius * radial_vector;
        
            delta_x = repmat(particles(1:2,:),1,1, n_active_gs) - ...
                      repmat(reshape(gs_position,2,1, n_active_gs),1,n_particles,1);
        
        
            distance = vecnorm(delta_x);
        
            radial_components = pagemtimes(reshape(radial_vector,1,2, n_active_gs), delta_x./distance);
            azimuthal_components = pagemtimes(reshape(az_vector,1,2, n_active_gs), delta_x./distance);
        
            beta = asin(azimuthal_components);
        
            

            measurements = [round(distance, obj.rounding(1,1)); round(beta, obj.rounding(2,1))];

            measurements = measurements + repmat(mvnrnd([0;0], diag([0.01, 0.000001]), n_particles)',1,1, n_active_gs);

            detection = radial_components >= 0;

        end


        function [measurements, detection] = get_3D_observations(obj, particles, time, active_gs)
        
            % particles: 6*N array
            % ground_station: 1*M array (long)
            % time (s)
        
            %   Output: 
            % Observations: 3*N*M, first row is distance, second is beta. Beta is the
            % angle moving ot the east
            %
            % detection: 1*N*M boolean, indicates whether particle N is seen by
            % ground station M

            n_particles = size(particles, 2);
        
            n_active_gs = sum(active_gs);

            i_earth = 0.4;            %Real: 0.4
            earth_radius = 6371;
        
            gs_lat = obj.gs_location(1,active_gs);
            gs_long = obj.gs_location(2,active_gs) + 2*pi/(24*3600) * time;
        
            radial_vector = [
                cos(gs_lat).*cos(gs_long);
                cos(i_earth).*cos(gs_lat).*sin(gs_long) - sin(i_earth).*sin(gs_lat);
                sin(i_earth).*cos(gs_lat).*sin(gs_long) + cos(i_earth).*sin(gs_lat)
            ];
        
            az_vector = [
                -sin(gs_long);
                cos(i_earth).*cos(gs_long);
                sin(i_earth).*cos(gs_long)
            ];
        
            north_vector = [
                -sin(gs_lat).*cos(gs_long);
                -cos(i_earth).*sin(gs_lat).*sin(gs_long) - sin(i_earth).*cos(gs_lat);
                -sin(i_earth).*sin(gs_lat).*sin(gs_long) + cos(i_earth).*cos(gs_lat)
            ];
        
            gs_position = earth_radius * radial_vector;
        

            delta_x = repmat(particles(1:3,:),1,1,n_active_gs) - ...
                      repmat(reshape(gs_position,3,1,n_active_gs),1,n_particles,1);
        
        
            distance = vecnorm(delta_x);
        
            radial_components = pagemtimes(reshape(radial_vector,1,3,n_active_gs), delta_x./distance);
            azimuthal_components = pagemtimes(reshape(az_vector,1,3,n_active_gs), delta_x./distance);
            north_components = pagemtimes(reshape(north_vector,1,3,n_active_gs), delta_x./distance); 
        
            alpha = asin(north_components);
        
            beta = asin(azimuthal_components);
        
            measurements = [round(distance, obj.rounding(1,1)); round(alpha, obj.rounding(2,1)) ;round(beta, obj.rounding(3,1))];
        
            detection = (radial_components >= 0);
        end

        function T = obtain_transform_matrix(obj, psi, theta, i)
            % Matrix that goes from from ground station centered frame to
            % Earth centered frame

            if nargin == 2
                % Two dimensions
                T = [cos(psi), -sin(psi); sin(psi), cos(psi)];
            else
                % Three dimensions
                T = [
                    cos(theta).*cos(psi), -sin(psi), -sin(theta).*cos(psi);
                    cos(i).*cos(theta).*sin(psi) - sin(i).*sin(theta), cos(i).*cos(psi), -cos(i).*sin(theta).*sin(psi) - sin(i).*cos(theta);
                    sin(i).*cos(theta).*sin(psi) + cos(i).*sin(theta), cos(i).*cos(psi), -sin(i).*sin(theta).*sin(psi) + cos(i).*cos(theta)
                ];
            end
        end

        function C = get_C(obj, x)
            C = [];
            if obj.three_dimensional
                T = obj.obtain_transform_matrix();
            end
        end

    end
end

