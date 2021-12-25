clc; clear; close all

addpath(genpath("."));

IC = [12791; 0; 0; 0; 5.58*1.1*cos(0.3); 5.58*1.1*sin(0.3)];

target = obtain_3D_motion(IC, [0; 20000]);

ground_stations = { 
    struct("location", [40.45*pi/180; -4.37*pi/180], "R", diag([1, 1/80000, 1/80000]));
    struct("location", [-35.4*pi/180; 148.98*pi/180], "R", diag([1, 1/80000, 1/80000]));
    struct("location", [35.34*pi/180; -116.87*pi/180], "R", diag([1, 1/80000, 1/80000]))
                };

dispersion_models = struct( ...
    "standard", eye(6)/1000, ...
    "no_sight", diag([1/10000000000, 1/10000000000, 1/10000000000, 0, 0, 0]) ...
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


figure

for i=1:20000

    my_filter = my_filter.step(1);

    % if mod(i,10) ~= 0
    %     continue
    % end

    time = i;

    gs_long =  2*pi/(24*3600) * time;

    
    X= my_filter.S(1,:);
    Y = my_filter.S(2,:);
    Z = my_filter.S(3,:);
    
    real_one = deval(target, time);

    scatter3(real_one(1), real_one(2),real_one(3), 20)
    hold on
    scatter3(X,Y,Z, 6, "filled")

    plot3(earth(1,:), earth(2,:), earth(3,:));
    scatter3(6371*cos(gs_long), 6371*sin(gs_long),0, 5, "g");
    plot3([6371, real_one(1)], ...
         [0, real_one(2)], ...
         [0, real_one(3)])
    hold off
    xlim([real_one(1)-200, real_one(1)+200])
    ylim([real_one(2)-200, real_one(2)+200])
    zlim([real_one(3)-200, real_one(3)+200])
    % xlim([-15000, 15000])
    % ylim([-15000, 15000])
    % zlim([-6000, 6000])
    daspect([1 1 1])
    drawnow
end
