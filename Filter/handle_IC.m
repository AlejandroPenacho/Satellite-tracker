function S = handle_IC(IC, n_particles)
    %HANDLE_IC Summary of this function goes here
    %   Detailed explanation goes here
    if IC.mode == "fixed_IC"
        S = repmat(IC.X,1,n_particles);
    
    elseif IC.mode == "normal"

        S = mvnrnd(IC.X', IC.cov, n_particles)';

    elseif IC.mode == "uniform"

        S_min = IC.X(:,1);
        S_max = IC.X(:,2);
        random_array = rand(size(IC.X,1), n_particles);
        S = S_min + random_array.* (S_max - S_min);

    elseif IC.mode == "3d-start"
        x_0 = IC.X;
        d = norm(x_0);
        v_std = sqrt(398600/d);
        angle = rand(1,n_particles)*2*pi;
        v_tan = v_std * (0.8 + 0.4*rand(1,n_particles));
        v_norm = (rand(1,n_particles)*1-2)*0.1*v_std;


    elseif IC.mode == "2d-start"
        x_0 = IC.X;
        d = norm(x_0);
        v_std = sqrt(398600/d);
        v = v_std * (1.2- 0.4*rand(1,n_particles));
        angle = 1 - 2*rand(1,n_particles);
        alpha = atan2(x_0(2), x_0(1));
        angle = angle + alpha + pi/2;
        S = [ ...
            repmat(x_0,1,n_particles);  ...
            v.*[cos(angle); sin(angle)] ...
            ];
    else

        error("Wrong IC mode");
    end
end

