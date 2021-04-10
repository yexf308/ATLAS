function [FPT_sim,exit_phi,t_final] = MFPT_butane2_sim(weighted_dd2, well_threshold, output, N_IC, phi_store, chart,chart_sim_parameter)

        d                   = chart_sim_parameter.d;
        dt_s                = chart_sim_parameter.dt_s;
        connectivity        = chart_sim_parameter.connectivity;
        
        connectivity_indices = cell(1,size(connectivity,1));
        for k = 1:size(connectivity,1)
         connectivity_indices{k} = find(connectivity(k,:));
        end
        
        sqrtdts              = sqrt(dt_s);
        [X_int_store]        = set_well(output,  N_IC, well_threshold);
        FPT_sim              = zeros(1, N_IC);
        exit_phi             = zeros(1, N_IC); % -1 means left point, +1 means right point
        
if true, rng('default'); rng(1); end
disp(['Use ATLAS Simulator with ', num2str(N_IC), ' initial points'])



tstart= tic;
 
parfor i = 1:N_IC 
            RAND_RESERVOID_MAX_SIZE = 10^8;
            RAND_RESERVOIR_SIZE     = 10^5;
            rand_reservoir          = randn(d,RAND_RESERVOIR_SIZE) * sqrtdts;
            rand_reservoir_index    = 1;
            
            X_curr                  = X_int_store(i,:);
            phi_curr                = atan2(X_curr(6),X_curr(4));
            [~,nearest]             = min(abs(phi_store-phi_curr));
            neigh                   = connectivity_indices{nearest};  
            
            step                    = 0;

            while phi_curr< well_threshold(2) && phi_curr > well_threshold(1)
                      [X_curr_proj , b,  ~  , H_hat, ~ , nearest, neigh ]     = weighted_dd2( X_curr, chart, neigh, nearest,connectivity_indices );
                      X_curr                                                  = X_curr_proj' + b * dt_s + H_hat * rand_reservoir(:,rand_reservoir_index);
                      X_curr                                                  = X_curr';  
                      phi_curr                                                = atan2(X_curr(6),X_curr(4)); 
                      rand_reservoir_index                                    = rand_reservoir_index+1;
                      step                                                    = step + 1;
                      if rand_reservoir_index > RAND_RESERVOIR_SIZE
                            if RAND_RESERVOIR_SIZE*2 < RAND_RESERVOID_MAX_SIZE
                                RAND_RESERVOIR_SIZE     = RAND_RESERVOIR_SIZE*2;
                            end
                            rand_reservoir          = randn(d,RAND_RESERVOIR_SIZE) * sqrtdts;
                            rand_reservoir_index    = 1;
                      end
            end
            
            exit_phi(i) = phi_curr;
            FPT_sim(i)  = step*dt_s;

end
t_final      = toc(tstart);
disp(['The time spent is ',num2str(t_final/3600), ' hours']);
disp(['MFPT in the ATLAS simulator is ', num2str(mean(FPT_sim)),'+-',num2str(std(FPT_sim)*1.96/sqrt(N_IC))])
        

end

