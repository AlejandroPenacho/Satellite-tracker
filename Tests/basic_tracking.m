clc; clear; close all

addpath(genpath("."));

target = obtain_2D_motion([6791; 0; 0; 7.66*1.1], [0; 20000]);


ground_stations = { 
    struct("location", 0, "R", diag([1, 1/8000]));
    struct("location", pi/4, "R", diag([1, 1/8000]))
                };

my_filter = ParticleFilter(10000, ...
                            [6791; 0; 0; 7.66*1.3], ...
                            false, ...
                            eye(4)/10, ...
                            ground_stations, ...
                            target);


earth = zeros(2,3);
for i=1:1000
    earth(:,i) = 6371*[cos(i*2*pi/1000); sin(2*pi*i/1000)];
end


figure

for i=1:20000

    my_filter = my_filter.step(1);

    time = i;

    gs_long =  2*pi/(24*3600) * time;

    
    X= my_filter.S(1,:);
    Y = my_filter.S(2,:);
    
    real_one = deval(target, time);

    scatter(real_one(1), real_one(2), 20)
    hold on
    scatter(X,Y, 6, "filled")

    plot(earth(1,:), earth(2,:));
    scatter(6371*cos(gs_long), 6371*sin(gs_long), 5, "g");
    plot([6371*cos(gs_long), real_one(1)], [6371*sin(gs_long), real_one(2)])
    plot([6371*cos(gs_long + pi/4), real_one(1)], [6371*sin(gs_long + pi/4), real_one(2)])
    hold off
    xlim([real_one(1)-100, real_one(1)+100])
    ylim([real_one(2)-100, real_one(2)+100])
    % xlim([0, 10000])
    % ylim([0, 10000])
    daspect([1 1 1])
    drawnow
end