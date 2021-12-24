

function observations = get_2D_observations(ground_station, particles, time)

    % particles: 4*N array
    % ground_station: 1*1 array (long)
    % time (s)

    % Output: 2*N, first row is alpha, second is beta. alpha is the angle
    % with respect to the radial vector, moving to the north. Beta is the
    % angle moving ot the east

    earth_radius = 6371;

    gs_long = ground_station(2) + 2*pi/(24*3600) * time;

    radial_vector = [
        cos(gs_long);
        sin(gs_long)
    ];

    az_vector = [
        -sin(gs_long);
        cos(gs_long)
    ];

    gs_position = earth_radius * radial_vector;

    delta_x = [
        particles(1,:) - gs_position(1);
        particles(2,:) - gs_position(2);
    ];

    distance = vecnorm(delta_x);

    beta = asin(az_vector'*delta_x./distance);

    observations = [distance; beta];
end