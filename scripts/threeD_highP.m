clc; clear; close all

addpath(genpath("."));

N = 20500;

X_0 =  [12791; 0; 0; 0; 5.58*1.1*cos(0.3); 5.58*1.1*sin(0.3)];

IC = struct("X", X_0, ...
            "mode", "fixed_IC");

target = obtain_3D_motion(X_0, [0; N*1.2], [true, true]);

distance_precision = 0.1; %km
angle_precision = 0.000029; %rad, equal to 6 sec (1 min is too low)


ground_stations = { 
    struct("location", [40.45*pi/180; -4.37*pi/180], "precision", [distance_precision; angle_precision; angle_precision]);
    struct("location", [-35.4*pi/180; 148.98*pi/180], "precision", [distance_precision; angle_precision; angle_precision]);
    struct("location", [35.34*pi/180; -116.87*pi/180], "precision", [distance_precision; angle_precision; angle_precision])
                };

dispersion_models = struct( ...
    "standard", eye(6)/100000, ...
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




% for i=1:5665
%     my_filter = my_filter.step(1);
%     waitbar(i/5500)
% end

real_state = zeros(6,N);
estimated_state = zeros(6,N);
covariance = zeros(6,N);


figure
for i=1:N

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

    xlabel("X (km)")
    ylabel("Y (km)")
    zlabel("Z (km)")

    real_one = deval(target, my_filter.time);

    real_state(:,i) = real_one;
    [mu, cov] = my_filter.get_estimation();
    estimated_state(:,i) = mu;
    covariance(:,i) = [cov(1,1); cov(2,2); cov(3,3); cov(1,2); cov(2,3); cov(1,3)];
end

time = 0:(N-1);
error = estimated_state-real_state;

distance_error = sqrt(sum(error(1:3,:).^2,1));

std_dev = sqrt(covariance);

%% Plot

figure
subplot(2,2,[1 2])
plot(time, distance_error, "LineWidth", 1.5)
grid minor
xlim([0 time(end)])
xlabel("Time (s)")
ylabel("Error (km)")
set(gca, "FontSize", 16)

subplot(2,2,3)
plot(time(1500:1600), distance_error(1500:1600), "LineWidth", 1.5)
grid minor
xlabel("Time (s)")
ylabel("Error (km)")
title("2 GS")
set(gca, "FontSize", 16)

subplot(2,2,4)
plot(time(1800:1900), distance_error(1800:1900), "LineWidth", 1.5)
grid minor
xlabel("Time (s)")
ylabel("Error (km)")
title("1 GS")
set(gca, "FontSize", 16)


figure
subplot(2,2,[1 2])
plot(time, std_dev(1,:), "LineWidth", 1.5)
hold on
plot(time, std_dev(2,:), "LineWidth", 1.5)
plot(time, std_dev(3,:), "LineWidth", 1.5)
xlim([0 time(end)])
hold off
grid minor
xlabel("Time (s)")
ylabel("Std. deviation (km)")
legend("\sigma_{1}", "\sigma_{2}", "\sigma_{3}")
set(gca, "FontSize", 16)

subplot(2,2,3)
plot(time(1500:1600), std_dev(1,1500:1600), "LineWidth", 1.5)
hold on
plot(time(1500:1600), std_dev(2,1500:1600), "LineWidth", 1.5)
plot(time(1500:1600), std_dev(3,1500:1600), "LineWidth", 1.5)
hold off
grid minor
xlabel("Time (s)")
ylabel("Std. deviation (km)")
title("2 GS")
set(gca, "FontSize", 16)

subplot(2,2,4)
plot(time(1800:1900), std_dev(1,1800:1900), "LineWidth", 1.5)
hold on
plot(time(1800:1900), std_dev(2,1800:1900), "LineWidth", 1.5)
plot(time(1800:1900), std_dev(3,1800:1900), "LineWidth", 1.5)
hold off
grid minor
xlabel("Time (s)")
ylabel("Std. deviation (km)")
title("1 GS")
set(gca, "FontSize", 16)