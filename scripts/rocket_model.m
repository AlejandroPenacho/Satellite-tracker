clc; clear; close all

addpath(genpath("."));

target_0 = [6791; 0; 0; 7.66*1.1];

N = 300;

target = obtain_rocket_motion([6791; 0; 0; 7.66*1.1], [0; 8000]);

plot_rocket(target);

 %%

dispersion_models_high_p = struct( ...
    "standard", eye(4)/100000, ...
    "no_sight", zeros(4), ...
    "recovery", diag([1, 1, 1/1000, 1/1000]), ...
    "first_contact", diag([1/10, 1/10, 1, 1]) ...
    );

dispersion_models_low_p = struct( ...
    "standard", eye(4)/10, ...
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



high_p_filter = ParticleFilter(10000, ...
                            IC, ...
                            false, ...
                            dispersion_models_high_p, ...
                            ground_stations, ...
                            target);

low_p_filter = ParticleFilter(10000, ...
                            IC, ...
                            false, ...
                            dispersion_models_low_p, ...
                            ground_stations, ...
                            target);

earth = zeros(2,3);
for i=1:1000
    earth(:,i) = 6371*[cos(i*2*pi/1000); sin(2*pi*i/1000)];
end


figure

high_p_filter.plot_state();

real_state = zeros(4,N);
estimated_state_hp = zeros(4,N);
covariance_hp = zeros(2,N);

estimated_state_lp = zeros(4,N);
covariance_lp = zeros(2,N);

mahala = zeros(1,N);

for i=1:N

    high_p_filter = high_p_filter.step(1);
    low_p_filter = low_p_filter.step(1);
    % fprintf("%f\n", my_filter.filter_state.mean_mahalanobis);

    mahala(i) = high_p_filter.filter_state.mean_mahalanobis;

    high_p_filter.plot_state();

    real_one = deval(target, high_p_filter.time);

    xlim([real_one(1)-100, real_one(1)+100])
    ylim([real_one(2)-100, real_one(2)+100])

    title("2-D particle filter")
    xlabel("X (km)")
    ylabel("Y (km)")
    % xlim([-10000, 10000])
    % ylim([-10000, 10000])
    drawnow

    real_state(:,i) = real_one;
    [mu_hp, cov_hp] = high_p_filter.get_estimation();
    estimated_state_hp(:,i) = mu_hp;
    covariance_hp(:,i) = eig(cov_hp(1:2,1:2));

    [mu_lp, cov_lp] = low_p_filter.get_estimation();
    estimated_state_lp(:,i) = mu_lp;
    covariance_lp(:,i) = eig(cov_lp(1:2,1:2));    

end

time = 0:(N-1);
error_hp = estimated_state_hp -real_state;
distance_error_hp = sqrt(error_hp(1,:).^2 + error_hp(2,:).^2);

error_lp = estimated_state_lp -real_state;
distance_error_lp = sqrt(error_lp(1,:).^2 + error_lp(2,:).^2);

%% Plot

figure
subplot(2,1,1)
plot(time, distance_error_hp, "LineWidth", 1.5)
hold on
plot(time, distance_error_lp, "LineWidth", 1.5)
hold off
grid minor
xlabel("Time (s)")
ylabel("Error (km)")
legend("Q = 10^{-6}", "Q=10^{-1}")
set(gca, "FontSize", 16)


subplot(2,1,2)
plot(time, sqrt(covariance_hp(1,:)), "LineWidth", 1.5, "color", [0, 0.4470, 0.7410])
hold on
plot(time, sqrt(covariance_hp(2,:)), "--", "LineWidth", 1.5, "color", [0, 0.4470, 0.7410])
plot(time, sqrt(covariance_lp(1,:)), "LineWidth", 1.5, "color", [0.8500, 0.3250, 0.0980])
plot(time, sqrt(covariance_lp(1,:)), "--", "LineWidth", 1.5, "color", [0.8500, 0.3250, 0.0980])
hold off
grid minor
xlabel("Time (s)")
ylabel("Std. deviation (km)")
set(gca, "FontSize", 16)
legend("\sigma_{1}", "\sigma_{2}")
ylim([0 5])


figure
plot(time, mahala, "LineWidth", 1.5)
grid minor



function plot_rocket(ode)

        % Earth
    earth = zeros(2,1001);
    for i=1:1001
        earth(:,i) = 6371*[cos(i*2*pi/1000); sin(2*pi*i/1000)];
    end

    plot(earth(1,:), earth(2,:), "color", [0.6863, 0.3690, 0], "LineWidth", 1.5)

    hold on

    time_1=0:300;
    time_2=300:700;

    X_1 = deval(ode, time_1)';
    X_2 = deval(ode, time_2)';

    h2 = plot(X_1(:,1), X_1(:,2), "Color", [0.8500, 0.3250, 0.0980], "LineWidth", 1.5);
    hold on
    plot(X_2(:,1), X_2(:,2), "--", "Color", [0.8500, 0.3250, 0.0980], "LineWidth", 1.5);
    h1 = scatter(6371, 0, "filled", "Color", [0, 0.4470, 0.7410]);
    hold off
    xlim([6200, 7000])
    ylim([-1000, 4000])
    xlabel("X (km)")
    ylabel("Y (km)")

    daspect([1 1 1])

    legend([h1, h2], ["Ground station", "Rocket traj."])

    set(gca, "FontSize", 15);

end