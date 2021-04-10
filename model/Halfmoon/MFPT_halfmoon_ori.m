function [FPT_ori, t_final] = MFPT_halfmoon_ori(well_threshold, output, phi_ori,N_IC,RHS_parameter,chart_sim_parameter, nonlin_trans_inv)
       drift                 = RHS_parameter.drift;
       diffusion             = RHS_parameter.diffusion;
       dt                    = RHS_parameter.dt;
       D                     = RHS_parameter.D;
       dt_s                  = chart_sim_parameter.dt_s;
       ratio                 = dt_s/dt;
       
       if well_threshold(2)>pi
           phi_ori ( phi_ori < (well_threshold(2)-2*pi) ) =  phi_ori ( phi_ori < (well_threshold(2)-2*pi) ) + 2*pi;
       end
       
       IC                    = output( phi_ori < well_threshold(2) & phi_ori > well_threshold(1),   :);
       [L_IC,~]              = size(IC);
       IC                    = IC(round(linspace(1, L_IC, N_IC)),:);
       
     
       FPT_ori               = zeros(1, N_IC);
       disp(['analysis in the well [', num2str(well_threshold(1)), ',',num2str(well_threshold(2)),']'])
       disp(['Use Original Simulator with ', num2str(N_IC), ' initial points'])
       
       sqrtdt                = sqrt(dt);
       tstart                = tic;
       
       if true, rng('default'); rng(1); end
       parfor i = 1: N_IC
           
            X           = IC(i,:)';
            Current     = nonlin_trans_inv(X);
            phi         = Current(1);
            step        = 0;
            while phi< well_threshold(2) && phi > well_threshold(1) || mod(step, ratio) ~=0
                drift_current         = drift(Current) * dt;
                diffusion_current     = diffusion(Current) * randn(D, 1) * sqrtdt;
                Current               = Current +  drift_current + diffusion_current;
                phi                   = Current(1);
                step = step+1;
            end
            FPT_ori(i)  = step*dt;
       end
t_final      = toc(tstart);
disp(['The time spent is ',num2str(t_final/3600), ' hours']); 
disp(['MFPT in the original simulator is ', num2str(mean(FPT_ori)),'+-',num2str(std(FPT_ori)*1.96/sqrt(N_IC))])
       
       
end

