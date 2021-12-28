function Z = obtain_3D_motion(X_0, time_interval, active_perturbations)


    % Active perturbations:
    % J2, J3

    mu_earth = 398600;

    Z = ode45(@(t,x) simple_3D_orbital_model(t,x, mu_earth, active_perturbations), time_interval, X_0);
    
    
end


function x_dot = simple_3D_orbital_model(~,x, mu_earth, active_perturbations)

    %
    %   x(1): x-position
    %   x(2): y-position
    %   x(3): z-position
    %   x(4): x-velocity
    %   x(5): y-velocity
    %   x(6): z-velocity

    J_2 = 1.75553e10;
    J_3 = -2.61913e11;
    
    r = (x(1)^2 + x(2)^2 + x(3)^2)^0.5;

    force_factor = mu_earth / r^3;
    
    x_dot = [
        x(4);
        x(5);
        x(6);
        - force_factor * x(1);
        - force_factor * x(2);
        -force_factor * x(3)
        ];
    
    if active_perturbations(1)
        x_dot = x_dot + [
            0; 0; 0;
            (J_2*x(1)/r^7)*(6*x(3)^2 - (3/2)*(x(1)^2+x(2)^2));
            (J_2*x(2)/r^7)*(6*x(3)^2 - (3/2)*(x(1)^2+x(2)^2));
            (J_2*x(3)/r^7)*(3*x(3)^2 - (9/2)*(x(1)^2+x(2)^2))
            ];
    end

    if active_perturbations(2)
        x_dot = x_dot + [
            0; 0; 0;
            (J_3*x(1)*x(3)/r^9) * (10*x(3)^2- (15/2)*(x(1)^2+x(2)^2));
            (J_3*x(2)*x(3)/r^9) * (10*x(3)^2- (15/2)*(x(1)^2+x(2)^2));
            (J_3/r^9) * (4*(x(3)^2)*(x(3)^2 - 3*(x(1)^2+x(2)^2)) + (3/2)*(x(1)^2+x(2)^2)^2)            
            ];
    end

end