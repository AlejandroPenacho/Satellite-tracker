clc; clear; close all

addpath(genpath("."));

X_0 =  [12791; 0; 0; 0; 5.58*1.1*cos(0.3); 5.58*1.1*sin(0.3)];

IC = struct("X", X_0, ...
            "mode", "fixed_IC");

target = obtain_3D_motion(X_0, [0; 300000], [true, true]);

distance_precision = 0.1; %km
angle_precision = 0.000029; %rad, equal to 6 sec (1 min is too low)


ground_stations = { 
    struct("location", [40.45*pi/180; -4.37*pi/180], "precision", [distance_precision; angle_precision; angle_precision]);
    struct("location", [-35.4*pi/180; 148.98*pi/180], "precision", [distance_precision; angle_precision; angle_precision]);
    struct("location", [35.34*pi/180; -116.87*pi/180], "precision", [distance_precision; angle_precision; angle_precision])
                };

dispersion_models = struct( ...
    "standard", eye(6)/1000, ...
    "no_sight", diag([0, 0, 0, 0, 0, 0]), ...
    "recovery", diag([20, 20, 20, 0.001, 0.001, 0.001]), ...
    "first_contact", diag([1/10, 1/10, 1/10, 1, 1, 1]) ...
    );

my_filter = ParticleFilter(10000, ...
                            IC, ...
                            true, ...
                            dispersion_models, ...
                            ground_stations, ...
                            target);


earth = zeros(3,1000);
for i=1:1000
    earth(:,i) = 6371*[cos(i*2*pi/1000); sin(2*pi*i/1000); 0];
end




for i=1:5665
    my_filter = my_filter.step(1);
    waitbar(i/5500)
end

figure
for i=1:20000

    my_filter = my_filter.step(1);

    my_filter.plot_state();

    X = deval(my_filter.target, my_filter.time);
    xlim([X(1) - 100, X(1) + 100])
    ylim([X(2) - 100, X(2) + 100])
    zlim([X(3) - 100, X(3) + 100])
    % xlim([-15000, 15000])
    % ylim([-15000, 15000])
    % zlim([-6000, 6000])
    drawnow


%     % if mod(i,10) ~= 0
%     %     continue
%     % end
% 
%     time = i;
% 
%     gs_long =  2*pi/(24*3600) * time;
% 
%     
%     X= my_filter.S(1,:);
%     Y = my_filter.S(2,:);
%     Z = my_filter.S(3,:);
%     
%     real_one = deval(target, time);
% 
%     scatter3(real_one(1), real_one(2),real_one(3), 20)
%     hold on
%     scatter3(X,Y,Z, 6, "filled")
% 
%     plot3(earth(1,:), earth(2,:), earth(3,:));
%     scatter3(6371*cos(gs_long), 6371*sin(gs_long),0, 5, "g");
%     plot3([6371, real_one(1)], ...
%          [0, real_one(2)], ...
%          [0, real_one(3)])
%     hold off
%     % xlim([real_one(1)-50, real_one(1)+50])
%     % ylim([real_one(2)-50, real_one(2)+50])
%     % zlim([real_one(3)-50, real_one(3)+50])
%     xlim([-15000, 15000])
%     ylim([-15000, 15000])
%     zlim([-6000, 6000])
%     daspect([1 1 1])
%     drawnow
end
