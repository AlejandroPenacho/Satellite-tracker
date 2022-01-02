classdef Updater
    % Represents the update operation of the filter, taking the set of
    % particles S_bar and the measurements from the ground stations. Then,
    % it uses them in order to obtain a new set of particles S.

    properties
        n_particles
        R
        R_inv
        three_dimensional
    end

    methods
        function obj = Updater(n_particles, three_dimensional)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here

            % if three_dimensional
            %     assert(size(R,1) == 3 && size(R,2) == 3)
            % else
            %     assert(size(R,1) == 2 && size(R,2) == 2)
            % end

            obj.n_particles = n_particles;

            obj.R = [];

            obj.R_inv = obj.R^(-1);
            obj.three_dimensional = three_dimensional;
        end

        function [S, filter_state] = update(obj, S_bar, real_obs, particle_obs, R, filter_state)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            if isempty(real_obs)
                S = S_bar;

            else
                weights = obj.obtain_weights(real_obs, particle_obs, R);
                filter_state.weight_variance = sum((weights- 1/obj.n_particles).^2)/obj.n_particles;
                S = obj.resample(S_bar, weights);
                
            end

        end

        function weights = obtain_weights(obj, real_obs, particle_obs, R)
           

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


            reshaped_real_obs = reshape(real_obs,[],1);

            reshaped_particle_obs = reshape(permute(particle_obs, [1,3,2]),[],obj.n_particles);

            Z = reshaped_particle_obs - repmat(reshaped_real_obs, 1, obj.n_particles);
            Z(2:end,:) = mod(Z(2:end,:) + pi, 2*pi) - pi;

            while true
                Z_mod = linsolve(R,Z);
    
                delta = sum(reshape(reshape(Z,1, []) .* reshape(Z_mod,1,[]), [], obj.n_particles) ,1);
    
                % delta = sum(reshape(reshape(Z,1, 2*obj.n_particles) .* reshape(Z_mod,1,2*obj.n_particles), 2, obj.n_particles) ,1);
        
                % TODO: Determinant should be updated!!!
    
                psi = (1/(2*pi*sqrt(det(obj.R))))*exp(-0.5*delta);

                if sum(psi > 10^(-3)) > 100
                    break
                else
                    R = 10*R;
                end
            end

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