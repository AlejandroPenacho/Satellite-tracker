clc; clear; close all

addpath(genpath("."));

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

my_filter = ParticleFilter(10000, ...
                            [6791; 0; 0; 7.66*1.3], ...
                            false, ...
                            dispersion_models, ...
                            ground_stations, ...
                            target);


earth = zeros(2,3);
for i=1:1000
    earth(:,i) = 6371*[cos(i*2*pi/1000); sin(2*pi*i/1000)];
end


figure

for i=1:40000

    my_filter = my_filter.step(1);

    my_filter.plot_state();

    real_one = deval(target, my_filter.time);

    xlim([real_one(1)-100, real_one(1)+100])
    ylim([real_one(2)-100, real_one(2)+100])

    % xlim([-10000, 10000])
    % ylim([-10000, 10000])
    drawnow

end