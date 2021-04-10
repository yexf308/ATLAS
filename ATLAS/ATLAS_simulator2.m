function [X, nearest_store, chart] = ATLAS_simulator2(weighted_dd2, chart_sim_parameter,RHS_parameter, simulator_par, chart)

X_int                    = chart_sim_parameter.X_int;
dt_s                     = chart_sim_parameter.dt_s;
D                        = chart_sim_parameter.D;
d                        = chart_sim_parameter.d;
Nstep                    = chart_sim_parameter.Nstep;
nearest                  = chart_sim_parameter.nearest;
connectivity             = chart_sim_parameter.connectivity;
gap                      = chart_sim_parameter.gap;
explore_threshold        = chart_sim_parameter.explore_threshold;

connectivity_indices     = cell(1,size(connectivity,1));
for k  = 1:size(connectivity,1)
        connectivity_indices{k} = find(connectivity(k,:));
end

T_max                     = RHS_parameter.T_max;
N                         = RHS_parameter.N;
dt                        = RHS_parameter.dt;
t0                        = RHS_parameter.t0;
chi_p                     = RHS_parameter.chi_p;
threshold                 = RHS_parameter.threshold;
modify                    = RHS_parameter.modify;
connectivity_threshold    = RHS_parameter.connectivity_threshold;


X                         = zeros( round(Nstep/gap), D); 
neigh                     = connectivity_indices{nearest};
X_curr                    = X_int;
index_store               = 1;
r_store                   = randn(d, Nstep) * sqrt(dt_s);
nearest_store             = zeros(1,Nstep/gap);

for j = 1:Nstep
     [X_curr_proj , b,  ~  , H_hat, T , nearest, neigh ]       = weighted_dd2( X_curr, chart, neigh, nearest,connectivity_indices );
     if sum(T < explore_threshold)==0 
         disp(['At step ',num2str(j),' add the landmark ',num2str(length(chart)+1), ' to the ATLAS'])
         Sim_parameter                = struct(                          ...
                               'T_max',        T_max,                    ...
                               'dt',           dt ,                      ...
                               'N',            N,                        ...
                               'X_int',        X_curr_proj               ...
                               );
       [data, Cov_store, Mean_store] = simulator_par(Sim_parameter);
       [ chart_add ]                 = Learning_Slow_Manifold(data, Cov_store, Mean_store, D, d,modify);
       [ chart, connectivity,  ~]    = landmark_add(chart, connectivity, chart_add,t0, chi_p, threshold(1), connectivity_threshold); % append newest chart to the end of chart lists    
       connectivity_indices          = cell(1,size(connectivity,1));
       for k  = 1:size(connectivity,1)
                 connectivity_indices{k} = find(connectivity(k,:));
       end
       nearest                       = length(chart);  % the added landmark is the last one in the chart
       neigh                         = connectivity_indices{nearest};      
       X_curr_proj                   = chart_add.X_int;
       b                             = chart_add.b;
       H_hat                         = chart_add.U * diag(chart_add.sigma);
       
     end  
     
     X_curr                                                   = X_curr_proj' + b * dt_s + H_hat * r_store(:,j);
     X_curr                                                   = X_curr';
     if mod(j,gap) == 0
         X( index_store,:)                             = X_curr_proj;
         nearest_store(index_store)                    = nearest;
         index_store                                   = index_store +1; 
     end

end

end

