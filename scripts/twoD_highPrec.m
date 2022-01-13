clc; clear; close all

addpath(genpath("."));

target_0 = [6791; 0; 0; 7.66*1.1];

N = 8000;

target = obtain_2D_motion([6791; 0; 0; 7.66*1.1], [0; 8000]);

dispersion_models = struct( ...
    "standard", eye(4)/100000, ...
    "no_sight", zeros(4), ...
    "recovery", diag([1, 1, 1/1000, 1/1000]), ...
    "first_contact", diag([1/10, 1/10, 1, 1]) ...
    );

distance_precision = 0.1; %km
angle_precision = 0.00029; %rad, equal to 1 sec

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

real_state = zeros(4,N);
estimated_state = zeros(4,N);
covariance = zeros(3,N);

mahala = zeros(1,N);

for i=1:N

    my_filter = my_filter.step(1);
    % fprintf("%f\n", my_filter.filter_state.mean_mahalanobis);

    mahala(i) = my_filter.filter_state.mean_mahalanobis;

    my_filter.plot_state();

    real_one = deval(target, my_filter.time);

    xlim([real_one(1)-1000, real_one(1)+1000])
    ylim([real_one(2)-1000, real_one(2)+1000])

    title("2-D particle filter")
    xlabel("X (km)")
    ylabel("Y (km)")
    % xlim([-10000, 10000])
    % ylim([-10000, 10000])
    drawnow

    real_state(:,i) = real_one;
    [mu, cov] = my_filter.get_estimation();
    estimated_state(:,i) = mu;
    covariance(:,i) = [cov(1,1); cov(2,2); cov(1,2)];

end

time = 0:(N-1);
error = estimated_state-real_state;
cov_comp = zeros(2,N);

distance_error = sqrt(error(1,:).^2 + error(2,:).^2);

for i=1:N
    cov_comp(:,i) = eig([covariance(1,i), covariance(3,i); covariance(3,i), covariance(2,i)]);
end

%% Plot

figure
subplot(2,1,1)
plot(time, distance_error, "LineWidth", 1.5)
grid minor
xlabel("Time (s)")
ylabel("Error (km)")
set(gca, "FontSize", 16)

axes("Position", [0.4 0.8 0.1 0.1])
box on
plot(time, distance_error, "LineWidth", 1.5)
xlim([2000, 2050])
ylim([0, 1])
set(gca, "FontSize", 16)
set(gca,'XTickLabel',[]);
grid on


subplot(2,1,2)
plot(time, sqrt(cov_comp(1,:)), "LineWidth", 1.5)
hold on
plot(time, sqrt(cov_comp(2,:)), "LineWidth", 1.5)
hold off
grid minor
xlabel("Time (s)")
ylabel("Std. deviation (km)")
set(gca, "FontSize", 16)
legend("\sigma_{1}", "\sigma_{2}")

axes("Position", [0.4 0.3 0.1 0.1])
box on
plot(time, sqrt(cov_comp(1,:)), "LineWidth", 1.5)
hold on
plot(time, sqrt(cov_comp(2,:)), "LineWidth", 1.5)
hold off
xlim([2000, 2050])
ylim([0, 1])
set(gca, "FontSize", 16)
set(gca,'XTickLabel',[]);
grid on


figure
plot(time, mahala, "LineWidth", 1.5)
grid minor