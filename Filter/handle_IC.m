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
    end
end

