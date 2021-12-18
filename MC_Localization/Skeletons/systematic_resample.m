% This function performs systematic re-sampling
% Inputs:   
%           S_bar(t):       4XM
% Outputs:
%           S(t):           4XM
function S = systematic_resample(S_bar)
	
    global M % number of particles 
    
    % YOUR IMPLEMENTATION
    
    S = zeros(4,M);
    

    if(isnan(S_bar(4,:)))
        S_bar(4,:) = (1/M)*ones(1,M);
    end

    CDF = cumsum(S_bar(4,:))/sum(S_bar(4,:));
    CDF(end) = 1;
    
    r = rand()/M;
    
    j = 1;
    for m=1:M
        while true
            if (r + (m-1)/M) <= CDF(j)
                S(:,m) = S_bar(:,j);
                break;
            else
                j = j+1;
            end
        end
    end
    
    
    
%     S = zeros(4,M);
%     
%     % YOUR IMPLEMENTATION
%      if(isnan(S_bar(4,:)))
%          S_bar(4,:) = (1/M)*ones(1,M);
%      end
%     
%     CDF = cumsum(S_bar(4,:));
%  
%     r = rand(1)*(1/M);
%     
%     for m = 1:M
%         i = find(CDF >= (r+(m-1)/M), 1, 'first');
%         
%         S(:,m) = S_bar(:,i);
%         S(4,m) = 1/M; % Recomputes weights
%     end

end