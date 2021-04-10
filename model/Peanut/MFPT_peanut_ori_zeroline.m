 function [FPT_ori, t_final] = MFPT_peanut_ori_zeroline(X_int_store,RHS_parameter,phi_zero, theta_zero)           
        drift               = RHS_parameter.drift;
        diffusion           = RHS_parameter.diffusion;
        dt                  = RHS_parameter.dt;
        t0                  = RHS_parameter.t0;
        ratio               = round(t0/dt);
        D                   = RHS_parameter.D;
        sqrtdt              = sqrt(dt);
        [N_IC, ~]           = size(X_int_store);
        FPT_ori             = zeros(1, N_IC);
        
         
if true, rng('default'); rng(1); end
        
disp(['Use Original Simulator with ', num2str(N_IC), ' initial points'])
        
tstart= tic;

parfor i = 1:N_IC 
        
         RAND_RESERVOID_MAX_SIZE = 10^8;
         RAND_RESERVOIR_SIZE     = 5*10^6;
         rand_reservoir          = randn(D,RAND_RESERVOIR_SIZE) * sqrtdt;
         rand_reservoir_index    = 1;
         X_curr                  = X_int_store(i,:)';
         step     = 0;
         ii       = 1;
         
         [index_int] = phitheta(X_curr, phi_zero, theta_zero); 
         
         while true
           
             [drift_curr]                                = drift( X_curr );
             [diffusion_curr]                            = diffusion( X_curr );
             drift_current                               = drift_curr * dt;
             diffusion_current                           = diffusion_curr * rand_reservoir(:,rand_reservoir_index) ;
             X_curr                                      = X_curr +  drift_current + diffusion_current;
             
             
    
               if mod(step, ratio)==0
                           [index_curr] = phitheta(X_curr, phi_zero, theta_zero);
                           if index_curr ~= index_int 
                                     break
                           end
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

t_final      = toc(tstart);
disp(['The time spent is ',num2str(t_final/3600), ' hours']);
disp(['MFPT in the original simulator is ', num2str(mean(FPT_ori)),'+-',num2str(std(FPT_ori)*1.96/sqrt(N_IC))])

        
end

