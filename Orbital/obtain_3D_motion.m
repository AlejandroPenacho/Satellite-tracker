function Z = obtain_3D_motion(X_0, time_interval)

    mu_earth = 398600;

    Z = ode45(@(t,x) simple_3D_orbital_model(t,x, mu_earth), time_interval, X_0);
    
    
end


function x_dot = simple_3D_orbital_model(~,x, mu_earth)

    %
    %   x(1): x-position
    %   x(2): y-position
    %   x(3): z-position
    %   x(4): x-velocity
    %   x(5): y-velocity
    %   x(6): z-velocity
    
    force_factor = mu_earth / (x(1)^2 + x(2)^2 + x(3)^2)^(3/2);
    
    x_dot = [
        x(4);
        x(5);
        x(6);
        - force_factor * x(1);
        - force_factor * x(2);
        -force_factor * x(3)
        ];
    

end