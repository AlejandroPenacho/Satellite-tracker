classdef Updater
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here


    % Ground station position is not actually calculated!
    properties
        n_particles
        R
        R_inv
        gs_location
        three_dimensional
    end

    methods
        function obj = Updater(n_particles, ground_stations, three_dimensional)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here

            % if three_dimensional
            %     assert(size(R,1) == 3 && size(R,2) == 3)
            % else
            %     assert(size(R,1) == 2 && size(R,2) == 2)
            % end

            obj.n_particles = n_particles;

            obj.R = [];
            obj.gs_location = [];
            for i=1:length(ground_stations)
                obj.R = blkdiag(obj.R, ground_stations{i}.R);
                obj.gs_location = [obj.gs_location, ground_stations{i}.location];
            end

            obj.R_inv = obj.R^(-1);
            obj.three_dimensional = three_dimensional;
        end

        function [S, update_data] = update(obj,S_bar, X, time, predict_data)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.three_dimensional
                [weights, detected] = obj.obtain_3D_weights(S_bar, X, time);
            else
                [weights, detected] = obj.obtain_2D_weights(S_bar, X, time);
            end

            if detected 
                S = obj.resample(S_bar, weights);
            else 
                S = S_bar;
            end

            update_data = struct("detected", detected);
        end

        function [weights, detected] = obtain_2D_weights(obj, S_bar, X, time)

            % n_gs = size(obj.gs_location, 2);
            [particle_obs, particle_detection] = get_2D_observations(obj.gs_location, S_bar, time);

            [real_obs, real_detection] = get_2D_observations(obj.gs_location, X, time);
           

            % One 3-dim for each ground station
            % 
            % Z = reshape(particle_obs(repmat(real_detection, 2, obj.n_particles,1)), 2 ,obj.n_particles,[]) - ...
            %     repmat(reshape(real_obs(repmat(real_detection, 2, 1,1)),2,1,[]),1, obj.n_particles,1);
            %
            % Only one 3-dim, additional ground stations in additional
            % rows. Non detections disappear
            %
            % Z = reshape(particle_obs(repmat(real_detection, 2, obj.n_particles,1)), [] ,obj.n_particles) - ...
            %     repmat(reshape(real_obs(repmat(real_detection, 2, 1,1)),[],1),1, obj.n_particles,1);

            % Only one 3-dim, measurements not detectable in real target
            % are set to 0, so they do not affect psi. R is the blkdiag of
            % the R of each ground station.

            detected = (sum(real_detection) > 0);

            Z = reshape(permute(particle_obs.*repmat(real_detection,2,obj.n_particles), [1,3,2]), [], obj.n_particles) - ...
                repmat(reshape(real_obs.*real_detection, [], 1),1, obj.n_particles,1);

            Z_mod = obj.R_inv * Z;

            delta = sum(reshape(reshape(Z,1, []) .* reshape(Z_mod,1,[]), [], obj.n_particles) ,1);

            % delta = sum(reshape(reshape(Z,1, 2*obj.n_particles) .* reshape(Z_mod,1,2*obj.n_particles), 2, obj.n_particles) ,1);
    
            % TODO: Detrminant should be updated!!!

            psi = (1/(2*pi*sqrt(det(obj.R))))*exp(-0.5*delta);

            weights = psi/sum(psi);

        end


        function [weights, detected] = obtain_3D_weights(obj, S_bar, X, time)

            % n_gs = size(obj.gs_location,2);

            [particle_obs, particle_detection] = get_3D_observations(obj.gs_location, S_bar, time);
            [real_obs, real_detection] = get_3D_observations(obj.gs_location, X, time);
           

            detected = (sum(real_detection) > 0);

            Z = reshape(permute(particle_obs.*repmat(real_detection,3,obj.n_particles), [1,3,2]), [], obj.n_particles) - ...
                repmat(reshape(real_obs.*real_detection, [], 1),1, obj.n_particles,1);

            Z_mod = obj.R_inv * Z;

            delta = sum(reshape(reshape(Z,1, []) .* reshape(Z_mod,1,[]), [], obj.n_particles) ,1);

            % TODO: Determinant should be updated!!!

            psi = (1/(2*pi*sqrt(det(obj.R))))*exp(-0.5*delta);

            weights = psi/sum(psi);

        end


        function S = resample(obj, S_bar, weights)

            if obj.three_dimensional 
                S = zeros(6, obj.n_particles);
            else
                S = zeros(4, obj.n_particles);
            end            
        
            if(isnan(weights))
                weights = (1/obj.n_particles)*ones(1,obj.n_particles);
            end
        
            CDF = cumsum(weights)/sum(weights);
            CDF(end) = 1;
            
            r = rand()/obj.n_particles;
            
            j = 1;
            for m=1:obj.n_particles
                while true
                    if (r + (m-1)/obj.n_particles) <= CDF(j)
                        S(:,m) = S_bar(:,j);
                        break;
                    else
                        j = j+1;
                    end
                end
            end

        end

    end
end