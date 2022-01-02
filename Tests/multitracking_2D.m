clc; clear; close all

addpath(genpath("."));

target_0 = [6791; 0; 0; 7.66];
target_1 = [-4500; 5362; -5.78; -4.85094];
target_2 = [4500; -5362; 5.78; 4.85094]; 


targets = {
    struct("detection_time", 0,"solved_ode", obtain_2D_motion(target_0, [0; 40000]));
    struct("detection_time", 0,"solved_ode", obtain_2D_motion(target_1, [0; 40000]));
    struct("detection_time", 0,"solved_ode", obtain_2D_motion(target_2, [0; 40000]))
    };

dispersion_models = struct( ...
    "standard", eye(4)/10, ...
    "no_sight", zeros(4), ...
    "recovery", diag([1, 1, 1/1000, 1/1000]), ...
    "first_contact", diag([1/10, 1/10, 1, 1]) ...
    );

distance_precision = 0.1; %km
angle_precision = 0.00029; %rad, equal to 1 min

ground_stations = { 
    struct("location", 0, "precision", [distance_precision; angle_precision]);
    struct("location", 2*pi/3, "precision", [distance_precision; angle_precision]);
    struct("location", 4*pi/3, "precision", [distance_precision; angle_precision])
    };


multitarget = MultitargetFilter(5000, false, dispersion_models, ground_stations, targets);

IC0 = struct("X", target_0, "cov", diag([5,5,0.1,0.1]), "mode", "normal");
IC1 = struct("X", target_1, "cov", diag([5,5,0.1,0.1]), "mode", "normal");
IC2 = struct("X", target_2, "cov", diag([5,5,0.1,0.1]), "mode", "normal");


% multitarget.current_filters = {
%     ParticleFilter(10000, IC0, false, dispersion_models, ground_stations, false);
%     ParticleFilter(10000, IC1, false, dispersion_models, ground_stations, false);
%     ParticleFilter(10000, IC2, false, dispersion_models, ground_stations, false)
% };

figure

for i=1:10000
    multitarget = multitarget.step(1);
    multitarget.plot();

    daspect([1 1 1])
    xlim([6000, 7000])
    ylim([-200, 5000])
    drawnow
end
