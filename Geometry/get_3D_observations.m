
function [observations, detection] = get_3D_observations(ground_station, particles, time)

    % particles: 6*N array
    % ground_station: 1*M array (long)
    % time (s)

    %   Output: 
    % Observations: 3*N*M, first row is distance, second is beta. Beta is the
    % angle moving ot the east
    %
    % detection: 1*N*M boolean, indicates whether particle N is seen by
    % ground station M

    n_gs = size(ground_station, 2);
    n_particles = size(particles, 2);

    i_earth = 0;            %Real: 0.4
    earth_radius = 6371;

    gs_lat = ground_station(1,:);
    gs_long = ground_station(2,:) + 2*pi/(24*3600) * time;

    radial_vector = [
        cos(gs_lat).*cos(gs_long);
        cos(i_earth).*cos(gs_lat).*sin(gs_long) - sin(i_earth).*sin(gs_lat);
        sin(i_earth).*cos(gs_lat).*sin(gs_long) + cos(i_earth).*sin(gs_lat)
    ];

    az_vector = [
        -sin(gs_long);
        cos(i_earth).*cos(gs_long);
        sin(i_earth).*cos(gs_long)
    ];

    north_vector = [
        -sin(gs_lat).*cos(gs_long);
        -cos(i_earth).*sin(gs_lat).*sin(gs_long) - sin(i_earth).*cos(gs_lat);
        -sin(i_earth).*sin(gs_lat).*sin(gs_long) + cos(i_earth).*cos(gs_lat)
    ];

    gs_position = earth_radius * radial_vector;

    delta_x_q = [
        particles(1,:) - gs_position(1);
        particles(2,:) - gs_position(2);
        particles(3,:) - gs_position(3)
    ];

    delta_x = repmat(particles(1:3,:),1,1,n_gs) - ...
              repmat(reshape(gs_position,3,1,n_gs),1,n_particles,1);


    distance = vecnorm(delta_x);

    radial_components = pagemtimes(reshape(radial_vector,1,3,n_gs), delta_x./distance);
    azimuthal_components = pagemtimes(reshape(az_vector,1,3,n_gs), delta_x./distance);
    north_components = pagemtimes(reshape(north_vector,1,3,n_gs), delta_x./distance); 

    alpha = asin(north_components);

    beta = asin(azimuthal_components);

    observations = [distance;alpha; beta];

    detection = (radial_components >= 0);
end