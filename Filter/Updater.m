classdef Updater
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here


    % Ground station position is not actually calculated!
    properties
        n_particles
        R
        three_dimensional
    end

    methods
        function obj = Updater(n_particles, R, three_dimensional)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here

            if three_dimensional
                assert(size(R,1) == 6 && size(R,2) == 6)
            else
                assert(size(R,1) == 4 && size(R,2) == 4)
            end

            obj.n_particles = n_particles;
            obj.R = R;
            obj.three_dimensional = three_dimensional;
        end

        function S = update(obj,S_bar, X, time)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.three_dimensional
                weights = obj.obtain_3D_weights(S_bar, X, time);
            else
                weights = obj.obtain_2D_weights(S_bar, X, time);
            end
            S = obj.resample(S_bar, weights);
        end

        function weights = obtain_2D_weights(obj, S_bar, X, time)

            ground_station = [0; 0];

            real_observations = get_2D_observations(ground_station, X, time);
           
            Z = get_2D_observations(ground_station, S_bar, time) - repmat(real_observations, 1, obj.n_particles);

            R_inv_d = 0.01;
            R_inv_b = 10000;

            psi = (1/(2*pi*sqrt(det(obj.R))))*exp(-0.5*(R_inv_d*Z(1,:).^2+R_inv_b*Z(2,:).^2));

            weights = psi/sum(psi);

        end


        function weights = obtain_3D_weights(obj, S_bar, X, time)

            ground_station = [0; 0];

            real_observations = get_3D_observations(ground_station, X, time);
           
            Z = get_3D_observations(ground_station, S_bar, time) - repmat(real_observations, 1, obj.n_particles);

            R_inv_d = 0.01;
            R_inv_b = 10000;

            psi = (1/(2*pi*sqrt(det(obj.R))))*exp(-0.5*(R_inv_d*Z(1,:).^2+R_inv_b*Z(2,:).^2+R_inv_b*Z(3,:).^2));

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