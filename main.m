clc; clear; close all

addpath(genpath("."));

target = obtain_3D_motion([6791; 0; 0; 0; 7.66*1.3*cos(0.3); 7.66*1.3*sin(0.3)], [0; 20000]);

ground_stations = { 
    struct("location", [0; 0], "R", diag([1, 1/8000, 1/8000]));
    struct("location", [0; pi/8], "R", diag([1, 1/8000, 1/8000]))
                };

dispersion_models = struct( ...
    "standard", eye(6)/10, ...
    "no_sight", diag([1/10000000000, 1/10000000000, 1/10000000000, 0, 0, 0]) ...
    );

my_filter = ParticleFilter(10000, ...
                            [6791; 0; 0; 0; 7.66*1.3*cos(0.3); 7.66*1.3*sin(0.3)], ...
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
    plot3([6371*cos(pi/8), real_one(1)], ...
         [6371*sin(pi/8), real_one(2)], ...
         [0, real_one(3)])
    hold off
    xlim([real_one(1)-200, real_one(1)+200])
    ylim([real_one(2)-200, real_one(2)+200])
    zlim([real_one(3)-200, real_one(3)+200])
    % xlim([5000, 7000])
    % ylim([0, 1000])
    % zlim([0, 1000])
    daspect([1 1 1])
    drawnow
end
