
function [FPT_sim, t_final] = MFPT_peanut_sim_general(weighted_dd2, well_threshold, X_int_store, nearest_store, chart,chart_sim_parameter,chart_angle)

   
d                   = chart_sim_parameter.d;
dt_s                = chart_sim_parameter.dt_s;
connectivity        = chart_sim_parameter.connectivity;
        
connectivity_indices = cell(1,size(connectivity,1));
for k = 1:size(connectivity,1)
       connectivity_indices{k} = find(connectivity(k,:));
end
        
sqrtdts              = sqrt(dt_s);
[N_IC, ~]            = size(X_int_store);
FPT_sim              = zeros(1, N_IC);
chart_angle          = chart_angle';
                
if true, rng('default'); rng(1); end
        
        
disp(['Use ATLAS Simulator with ', num2str(N_IC), ' initial points'])
tstart= tic; 
parfor i = 1:N_IC 
            RAND_RESERVOID_MAX_SIZE = 10^8;
            RAND_RESERVOIR_SIZE     = 10^5;
            rand_reservoir          = randn(d,RAND_RESERVOIR_SIZE) * sqrtdts;
            rand_reservoir_index    = 1;
            
            X_curr                  = X_int_store(i,:);
            step                    = 0;
            nearest                 = nearest_store(i);
            neigh                   = connectivity_indices{nearest};  
            while true
               [X_curr_proj , b,  ~  , H_hat, ~ , nearest, neigh ]       = weighted_dd2( X_curr, chart, neigh, nearest,connectivity_indices );  
                X_curr                                      = X_curr_proj' + b * dt_s + H_hat * rand_reservoir(:,rand_reservoir_index);
                X_curr                                      = X_curr';  

                theta_curr                  = atan2( sqrt( X_curr(1)^2 + X_curr(2)^2 ), X_curr(3) );
                phi_curr                    = mod(atan2(X_curr(2), X_curr(1)),2*pi);
                angle                       = [phi_curr, theta_curr];
                nearest_test                =  dsearchn(chart_angle,angle);    
                if ~ismember(nearest_test, well_threshold)
                         break
                end

                 rand_reservoir_index                        = rand_reservoir_index+1;
                 step                                        = step + 1;

                 if rand_reservoir_index > RAND_RESERVOIR_SIZE
                                if RAND_RESERVOIR_SIZE*2 < RAND_RESERVOID_MAX_SIZE
                                    RAND_RESERVOIR_SIZE     = RAND_RESERVOIR_SIZE*2;
                                end
                                rand_reservoir          = randn(d,RAND_RESERVOIR_SIZE) * sqrtdts;
                                rand_reservoir_index    = 1;
                 end
             end
       FPT_sim(i) = step*dt_s;

end
t_final      = toc(tstart);
disp(['The time spent is ',num2str(t_final/3600), ' hours']);
disp(['MFPT in the original simulator is ', num2str(mean(FPT_sim)),'+-',num2str(std(FPT_sim)*1.96/sqrt(N_IC))])
end

