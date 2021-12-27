classdef ParticleFilter
    %PARTICLEFILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        target
        S
        n_particles
        time
        filter_state
        three_dimensional
        observation_system
        predictor
        updater
    end
    
    methods
        function obj = ParticleFilter(n_particles, IC, three_dimensional, dispersion_models, ground_stations, target)
            %PARTICLEFILTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.n_particles = n_particles;
            obj.three_dimensional = three_dimensional;

            obj.S = repmat(IC,1,n_particles);

            obj.time = 0;

            obj.target = target;

            obj.predictor = Predictor(n_particles, dispersion_models, three_dimensional);
            obj.observation_system = ObservationSystem(three_dimensional, ground_stations);
            obj.updater = Updater(n_particles, three_dimensional);

            obj.filter_state = struct("detected", true);
        end
        
        function obj = step(obj, delta_t)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            obj.time = obj.time + delta_t;
            X = deval(obj.target, obj.time);
            

            [real_obs, real_detection] = obj.observation_system.get_measurements(X, obj.time);

            obj.filter_state.active_gs = reshape(real_detection,1,[]) == 1;

            real_obs = real_obs(:,:,real_detection==1);

            [S_bar, obj.filter_state] = obj.predictor.predict(obj.S, delta_t, obj.filter_state);
            
            

            [particle_obs, ~] = obj.observation_system.get_measurements(S_bar, obj.time, obj.filter_state.active_gs);
           
            R = obj.observation_system.get_R(obj.filter_state.active_gs);

            [obj.S, obj.filter_state] = obj.updater.update(S_bar, real_obs, particle_obs, R, obj.filter_state);
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

