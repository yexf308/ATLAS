       
        disp('Starting exploration learning stage')

        % set to exploration mode
        mode                    = 1;
        weighted_dd2            = @(X0, chart,neigh,nearest, connectivity_indices) weighted_drift_diffusion2( X0, chart, neigh, nearest, connectivity_indices, t0, chi_p, D,d, threshold, option,mode );


        %% Use the exploration mode to further Learn 
        explore_round                    = 1;
        explore_round_max                = 100; 
        K                                = length(chart);
        random_start;     
        while bin_N~=1 && explore_round<explore_round_max
              disp(['Round ',num2str(explore_round),' of exploration'])
              [~, ~, chart]                            = ATLAS_simulator2(weighted_dd2, chart_sim_parameter,RHS_parameter, simulator_par, chart);        
             [chart, ~, ~, ~, ~]                       = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
              disp(['No. of landmarks after exploration is ', num2str(length(chart)) ])

             if relearn_option == 1
                  K_next                                   = length(chart);
                  index_to_learn                           = K+1: K_next;
                 [ chart ]                                 = relearn_chart(chart, relearn_parameter, index_to_learn);
                 [chart, connectivity, P, bin_N, bins]     = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
             end
             K                                         = length(chart);
             random_start;    
             explore_round                             = explore_round + 1;
        end

        save( chart_fileName ,'chart','connectivity','P')
        
        disp('exploration learning stage is completed')