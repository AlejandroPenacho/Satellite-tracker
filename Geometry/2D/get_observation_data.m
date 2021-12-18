function [observation] = get_observation_data(station_angle, satellite_state_array)

    % station_angle: angle of the station, in radians
    % satellite_state_array: 4*N array for N satellites

    station_position = 6371 * [cos(station_angle); sin(station_angle)];
    
    distances = sqrt((station_position(1) - satellite_state_array(1,:))^2 + ...
                    (station_position(2) - satellite_state_array(2,:))^2);
                
    abs_angles = atan2(satellite_state_array(2,:) - station_position(2), ...
                   satellite_state_array(1,:) - station_position(1));
               
    real_angles = abs_angles - station_angle;
                
end

function [observation] = get_observation_data_2D(station_angle, satellite_state_array)

    % station_angle: angle of the station, in radians
    % satellite_state_array: 4*N array for N satellites

    station_position = 6371 * [cos(station_angle); sin(station_angle)];
    
    distances = sqrt((station_position(1) - satellite_state_array(1,:))^2 + ...
                    (station_position(2) - satellite_state_array(2,:))^2);
                
    abs_angles = atan2(satellite_state_array(2,:) - station_position(2), ...
                   satellite_state_array(1,:) - station_position(1));
               
    real_angles = abs_angles - station_angle;
                
end

function [observation] = get_observation_data_3D(station_state, satellite_state_array)

    % station_angle: angle of the station, in radians
    % satellite_state_array: 6*N array for N satellites
    
    station_lat = station_state(1);
    station_long = station_state(2);

    station_position = 6371 * [
        cos(station_lat)*cos(station_long);
        cos(station_lat)*cos(station_long); 
        sin(station_lat)];
                
    delta_x = [
        station_position(1) - satellite_state_array(1,:);
        station_position(2) - satellite_state_array(2,:);
        station_position(3) - satellite_state_array(3,:)
        ];
    
    distances = sqrt(delta_x(1,:).^2 + delta_x(2,:).^2 + delta_x(3,:).^2);
               
    T = [
                
end
