function [FPT_ori,exit_phi,t_final ] = MFPT_butane2_ori(well_threshold, output, N_IC, RHS_parameter,chart_sim_parameter)
       drift               = RHS_parameter.drift;
       diffusion           = RHS_parameter.diffusion;
       dt                  = RHS_parameter.dt;
       dt_s                = chart_sim_parameter.dt_s;
       ratio               = round(dt_s/dt);
       D                   = RHS_parameter.D;
       sqrtdt              = sqrt(dt);
       [X_int_store]       = set_well(output,  N_IC, well_threshold);
       FPT_ori             = zeros(1, N_IC);
       exit_phi            = zeros(1, N_IC); % -1 means left point, +1 means right point
       
      tstart= tic;
 if true, rng('default'); rng(1); end
 parfor i = 1:N_IC 
        
         RAND_RESERVOID_MAX_SIZE = 10^8;
         RAND_RESERVOIR_SIZE     = 5*10^6;
         rand_reservoir          = randn(D,RAND_RESERVOIR_SIZE) * sqrtdt;
         rand_reservoir_index    = 1;
         X_curr                  = X_int_store(i,:)';
         phi_curr = atan2(X_curr(6),X_curr(4));
         step     = 0;
         
         
         while phi_curr< well_threshold(2) && phi_curr > well_threshold(1) || mod(step, ratio) ~=0
              drift_current              = drift(X_curr) * dt;
              diffusion_current          = diffusion(X_curr) * rand_reservoir(:,rand_reservoir_index) ;
              X_curr                     = X_curr +  drift_current + diffusion_current;
              phi_curr                   = atan2(X_curr(6),X_curr(4));
              step                       = step + 1;
              rand_reservoir_index       = rand_reservoir_index+1;
          
             
              if rand_reservoir_index > RAND_RESERVOIR_SIZE
                        if RAND_RESERVOIR_SIZE*2 < RAND_RESERVOID_MAX_SIZE
                                RAND_RESERVOIR_SIZE     = RAND_RESERVOIR_SIZE*2;
                        end
                            rand_reservoir          = randn(D,RAND_RESERVOIR_SIZE) * sqrtdt ;
                            rand_reservoir_index    = 1;
              end
              
         end
        exit_phi(i) = phi_curr;
        FPT_ori(i)  = step*dt;
 end

t_final      = toc(tstart);
disp(['The time spent is ',num2str(t_final/3600), ' hours']);

disp(['MFPT in the original simulator is ', num2str(mean(FPT_ori)),'+-',num2str(std(FPT_ori)*1.96/sqrt(N_IC))])

       
       
end

