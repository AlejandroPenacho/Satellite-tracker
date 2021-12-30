clc; clear; close all

addpath(genpath("."));

target_0 = [6791; 0; 0; 7.66*1.1];

target = obtain_2D_motion([6791; 0; 0; 7.66*1.1], [0; 40000]);

dispersion_models = struct( ...
    "standard", eye(4)/10, ...
    "no_sight", zeros(4), ...
    "recovery", diag([1/10, 1/10, 1/10, 1/10]) ...
    );

distance_precision = 0.1; %km
angle_precision = 0.00029; %rad, equal to 1 min

ground_stations = { 
    struct("location", 0, "precision", [distance_precision; angle_precision]);
    struct("location", 2*pi/3, "precision", [distance_precision; angle_precision]);
    struct("location", 4*pi/3, "precision", [distance_precision; angle_precision])
    };

% IC = struct("X", [6791; 0; 0; 7.66*1.3], ...
%             "mode", "normal", ...
%             "cov", diag([5000, 5000, 500, 500]));

IC = struct("X", [6400, 8000; -1000, 1000; -100, 100; 7.66*0.8, 7.66*1.2], ...
            "mode", "uniform");



my_filter = ParticleFilter(10000, ...
                            IC, ...
                            false, ...
                            dispersion_models, ...
                            ground_stations, ...
                            target);


earth = zeros(2,3);
for i=1:1000
    earth(:,i) = 6371*[cos(i*2*pi/1000); sin(2*pi*i/1000)];
end


figure

my_filter.plot_state();

for i=1:40000

    my_filter = my_filter.step(1);
    fprintf("%f\n", my_filter.filter_state.weight_variance);

    my_filter.plot_state();

    real_one = deval(target, my_filter.time);

    xlim([real_one(1)-1000, real_one(1)+1000])
    ylim([real_one(2)-1000, real_one(2)+1000])

    % xlim([-10000, 10000])
    % ylim([-10000, 10000])
    drawnow

end