function [TranM] = MSM(chart, connectivity, MSM_parameter, weighted_dd2, d)
disp('Starting MSM part')

% This function gives the transition matrix M = exp(P*step*t0) 
step     = MSM_parameter.step;
dt_s     = MSM_parameter.dt_s;
N        = MSM_parameter.N_state; 
K        = length(chart);
TranM    = zeros(K,K);

connectivity_indices = cell(1,size(connectivity,1));
for k = 1:size(connectivity,1)
    connectivity_indices{k} = find(connectivity(k,:));
end



parfor i = 1:K
    TranM_i       = zeros(1,K);
    X_int         = chart{i}.X_int;
    nearest_store = zeros(1,N);

    for j = 1:N
        X            = X_int;   
        nearest      = i;
        neigh        = connectivity_indices{nearest};
        for k = 1:step
             [X_proj , b,  ~  , H_hat, ~ , nearest, neigh ]       = weighted_dd2( X, chart, neigh, nearest,connectivity_indices ); 
              X                                                   = X_proj' + b * dt_s + H_hat * randn(d,1) * sqrt(dt_s);
              X                                                   = X';
        end
             [~ , ~,  ~  , ~, ~ , nearest, ~ ]       = weighted_dd2( X, chart, neigh, nearest,connectivity_indices );         
             nearest_store(j)                        = nearest;
    end
    a                       = unique(nearest_store);
    TranM_i(a)              = histc(nearest_store,a)./N;
    TranM(i,:)              = TranM_i;
        
end

disp('MSM stage is completed')
end
