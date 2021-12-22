

function observations = get_observations(ground_station, particles, time)

    % particles: 6*N array
    % ground_station: 2*1 array (lat, long)
    % time (s)

    % Output: 2*N, first row is alpha, second is beta. alpha is the angle
    % with respect to the radial vector, moving to the north. Beta is the
    % angle moving ot the east

    i_earth = 0;            %Real: 0.4
    earth_radius = 6371;

    gs_lat = ground_station(1);
    gs_long = ground_station(2) + 2*pi/(24*3600) * time;

    radial_vector = [
        cos(gs_lat)*cos(gs_long);
        cos(i_earth)*cos(gs_lat)*sin(gs_long) - sin(i_earth)*sin(gs_lat);
        sin(i_earth)*cos(gs_lat)*sin(gs_long) + cos(i_earth)*sin(gs_lat)
    ];

    az_vector = [
        -sin(gs_long);
        cos(i_earth)*cos(gs_long);
        sin(i_earth)*cos(gs_long)
    ];

    north_vector = [
        -sin(gs_lat)*cos(gs_long);
        -cos(i_earth)*sin(gs_lat)*sin(gs_long) - sin(i_earth)*cos(gs_lat);
        -sin(i_earth)*sin(gs_lat)*sin(gs_long) + cos(i_earth)*cos(gs_long)
    ];

    gs_position = earth_radius * radial_vector;

    delta_x = [
        particles(1,:) - gs_position(1);
        particles(2,:) - gs_position(2);
        particles(3,:) - gs_position(3)
    ];

    distance = vecnorm(delta_x);

    beta = asin(az_vector'*delta_x./distance);

    alpha = asin(north_vector'*delta_x./distance);

    observations = [distance;alpha; beta];
end