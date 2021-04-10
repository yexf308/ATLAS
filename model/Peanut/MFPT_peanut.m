function [FPT_ori,FPT_sim] = MFPT_peanut(weighted_dd, well_threshold, X_int_store, nearest_store, chart,RHS_parameter,chart_sim_parameter)
        driftdiffusion      = RHS_parameter.driftdiffusion;
        diffusionM          = RHS_parameter.diffusionM;
        dt                  = RHS_parameter.dt;
        D                   = RHS_parameter.D;
        d                   = RHS_parameter.d;
        dt_s                = chart_sim_parameter.dt_s;
        connectivity        = chart_sim_parameter.connectivity;
        
        connectivity_indices = cell(1,size(connectivity,1));
        for k = 1:size(connectivity,1)
         connectivity_indices{k} = find(connectivity(k,:));
        end
        
        sqrtdt      = sqrt(dt);
        sqrtdts     = sqrt(dt_s);
        [N_IC, ~]            = size(X_int_store);
        FPT_ori              = zeros(1, N_IC);
        FPT_sim              = zeros(1, N_IC);
        
        K = length(chart);
        chart_angle       = zeros(d, K);

        for i = 1:K
             X_int                 = chart{i}.X_int';
             theta_X_int           = atan2( sqrt( X_int(1)^2 + X_int(2)^2 ), X_int(3) );
             phi_X_int             = mod(atan2(X_int(2), X_int(1)),2*pi);
             chart_angle(:,i)      = [phi_X_int;  theta_X_int ];
        end
        
        
        
        
        
disp(['Use Original Simulator with ', num2str(N_IC), ' initial points'])
        
tic
parfor i = 1:N_IC 
        
         RAND_RESERVOID_MAX_SIZE = 10^8;
         RAND_RESERVOIR_SIZE     = 5*10^6;
         rand_reservoir          = randn(D,RAND_RESERVOIR_SIZE) * sqrtdt;
         rand_reservoir_index    = 1;
        
  
         X_curr                  = X_int_store(i,:)';
         step     = 0;
         ii       = 1;
        
        while ii
           
             [drift,diffusion]                           = driftdiffusion( X_curr );
             drift_current                               = drift * dt;
             diffusion_current                           = diffusion * rand_reservoir(:,rand_reservoir_index) ;
             X_bar                                       = X_curr +  drift_current + diffusion_current;
             
             diffusion_bar                               = diffusionM(X_bar) * rand_reservoir(:,rand_reservoir_index);
             X_curr                                      = X_curr + drift_current + (diffusion_current + diffusion_bar)/2;
             
             % check current point is in the well or not
             theta_X_curr           = atan2( sqrt( X_curr(1)^2 + X_curr(2)^2 ), X_curr(3) );
             phi_X_curr             = mod(atan2(X_curr(2), X_curr(1)),2*pi);
             angle                  = [phi_X_curr, theta_X_curr];
             nearest                =  dsearchn(chart_angle',angle);            
             if ismember(nearest, well_threshold)
                          ii = 1;
             else
                          ii = 0;
             end
              
             
             rand_reservoir_index                        = rand_reservoir_index+1;
             step                                        = step + 1;
             
              if rand_reservoir_index > RAND_RESERVOIR_SIZE
                            if RAND_RESERVOIR_SIZE*2 < RAND_RESERVOID_MAX_SIZE
                                RAND_RESERVOIR_SIZE     = RAND_RESERVOIR_SIZE*2;
                            end
                            rand_reservoir          = randn(D,RAND_RESERVOIR_SIZE) * sqrtdt ;
                            rand_reservoir_index    = 1;
              end
              
        end
        FPT_ori(i) = step*dt;
end   
toc

disp(['MFPT in the original simulator is ', num2str(mean(FPT_ori)), ' and standard deiviation is ', num2str(std(FPT_ori))])


%% ATLAS simulator

if true, rng('default'); rng(1); end
disp(['Use ATLAS Simulator with ', num2str(N_IC), ' initial points'])


tic 
parfor i = 1:N_IC 
            RAND_RESERVOID_MAX_SIZE = 10^8;
            RAND_RESERVOIR_SIZE     = 10^5;
            rand_reservoir          = randn(d,RAND_RESERVOIR_SIZE) * sqrtdts;
            rand_reservoir_index    = 1;
            
            X_curr                  = X_int_store(i,:);
            
            step                    = 0;
            nearest                 = nearest_store(i);
            neigh                   = connectivity_indices{nearest};  
            ii                      = 1;
            while ii
                      [X_curr_proj, b,  ~  , H_hat, T, ~ ]        = weighted_dd( X_curr, chart, neigh );
                      [~, T_index]                                = min(T);
                      nearest                                     = neigh(T_index);
                      neigh                                       = connectivity_indices{nearest};
                      X_curr                                      = X_curr_proj' + b * dt_s + H_hat * rand_reservoir(:,rand_reservoir_index);
                      X_curr                                      = X_curr';  
                      
             
                      if ismember(nearest, well_threshold)
                          ii = 1;
                      else
                          ii = 0;
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
toc

disp(['MFPT in the ATLAS simulator is ', num2str(mean(FPT_sim)), ' and standard deiviation is ', num2str(std(FPT_sim))])

        
end

