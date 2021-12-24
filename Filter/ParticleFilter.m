classdef ParticleFilter
    %PARTICLEFILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        target
        S
        n_particles
        time
        predictor
        updater
    end
    
    methods
        function obj = ParticleFilter(n_particles, IC, three_dimensional, Q, ground_stations, target)
            %PARTICLEFILTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.n_particles = n_particles;

            obj.S = repmat(IC,1,n_particles);

            obj.time = 0;

            obj.target = target;

            obj.predictor = Predictor(n_particles, Q, three_dimensional);
            obj.updater = Updater(n_particles, ground_stations, three_dimensional);
        end
        
        function obj = step(obj, delta_t)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            S_bar = obj.predictor.predict(obj.S, delta_t);

            obj.time = obj.time + delta_t;

            X = deval(obj.target, obj.time);

            obj.S = obj.updater.update(S_bar, X, obj.time);
        end
    end
end

