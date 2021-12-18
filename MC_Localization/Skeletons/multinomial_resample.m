% This function performs multinomial re-sampling
% Inputs:   
%           S_bar(t):       4XM
% Outputs:
%           S(t):           4XM
function S = multinomial_resample(S_bar)

    global M % number of particles
    
    % YOUR IMPLEMENTATION
    
    S = zeros(4,M);
        
    CDF = cumsum(S_bar(4,:))/sum(S_bar(4,:));
    CDF(end) = 1;

    for n=1:M
        r = rand();
        for j=1:M
            if r<=CDF(j)
                S(:,n) = S_bar(:,j);
                break
            end
        end
    end
    
end
