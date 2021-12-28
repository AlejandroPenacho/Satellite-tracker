function out = obtain_2D_motion(X_0, time_interval, mode)

    
    mu_earth = 398600;


    out = ode45(@(t,x) simple_2D_orbital_model(t,x, mu_earth), time_interval, X_0);

end


function x_dot = simple_2D_orbital_model(~,x, mu_earth)

    %
    %   x(1): x-position
    %   x(2): y-position
    %   x(3): x-velocity
    %   x(4): y-velocity
    
    force_factor = mu_earth / (x(1)^2+x(2)^2)^(3/2);
    
    x_dot = [
        x(3);
        x(4);
        - force_factor * x(1);
        - force_factor * x(2)
        ];
    

end
