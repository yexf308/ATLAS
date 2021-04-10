function [FPT_sim, t_final] = MFPT_halfmoon_sim(weighted_dd2, well_threshold, output, phi_ori, N_IC, phi_store, chart,chart_sim_parameter)
        d                    = chart_sim_parameter.d;  
        dt_s                 = chart_sim_parameter.dt_s;
        connectivity         = chart_sim_parameter.connectivity;
        t                    = chart_sim_parameter.t;
       
        connectivity_indices = cell(1,size(connectivity,1));
       for k = 1:size(connectivity,1)
         connectivity_indices{k} = find(connectivity(k,:));
       end
       
       
       if well_threshold(2)>pi
           phi_ori ( phi_ori < (well_threshold(2)-2*pi) )  =  phi_ori ( phi_ori < (well_threshold(2)-2*pi) ) + 2*pi;
           phi_store(phi_store < (well_threshold(2)-2*pi) )= phi_store(phi_store < (well_threshold(2)-2*pi) ) + 2*pi;
       end
       
       IC                   = output( phi_ori < well_threshold(2) & phi_ori > well_threshold(1),   :);
       [L_IC,~]             = size(IC);
       IC                   = IC(round(linspace(1, L_IC, N_IC)),:);
       phi_IC               = phi_ori(phi_ori < well_threshold(2) & phi_ori > well_threshold(1));
       phi_IC               = phi_IC(round(linspace(1, L_IC, N_IC)));
       FPT_sim              = zeros(1, N_IC);
       
             
       disp(['Use ATLAS Simulator with ', num2str(N_IC), ' initial points'])
       tstart               = tic;
       sqrtdts              = sqrt(dt_s);
       if true, rng('default'); rng(1); end
       
        parfor i = 1:N_IC
            [~,j]                   = min(abs(phi_store-phi_IC(i)));
            X                       = IC(i,:);
            phi                     = atan2(X(2),X(1)) - sqrt(X(1)^2+X(2)^2)-t;
            step                    = 0;
            nearest                 = j;
            neigh                   = connectivity_indices{nearest};   
            while phi < well_threshold(2) && phi > well_threshold(1)
                
               [X_proj , b,  ~  , H_hat, ~ , nearest, neigh ]       = weighted_dd2( X, chart, neigh, nearest,connectivity_indices );
                X                                                   = X_proj' + b * dt_s + H_hat * randn(d,1) * sqrtdts;
                X                                                   = X';
                phi                                                 = atan2(X(2),X(1))- sqrt(X(1)^2+X(2)^2)-t;
                if  phi<-pi/2 && well_threshold(2)>pi 
                     phi            = phi + 2*pi;
                end
                step                                                = step + 1;
            end
            FPT_sim(i)              = step*dt_s;

            
        end
       t_final                      = toc(tstart);
       disp(['The time spent is ',num2str(t_final/3600), ' hours']); 

disp(['MFPT in the ATLAS simulator is ', num2str(mean(FPT_sim)),'+-',num2str(std(FPT_sim)*1.96/sqrt(N_IC))])
       
       
       
       
end

