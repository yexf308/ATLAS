        disp('Starting exploration learning stage')

        % set to exploration mode
        mode                    = 1;
        weighted_dd2            = @(X0, chart,neigh,nearest, connectivity_indices) weighted_drift_diffusion2( X0, chart, neigh,nearest, connectivity_indices, params.RHS.t0, params.RHS.chi_p, params.RHS.D,params.RHS.d, params.RHS.threshold, params.learning.option,mode );


        %% Use the exploration mode to further Learn
        explore_round                    = 1;
        explore_round_max                = 100;
        K                                = length(chart);
        random_start;
        while bin_N~=1 && explore_round<explore_round_max
              disp(['Round ',num2str(explore_round),' of exploration'])
              [~, ~, chart]                            = ATLAS_simulator2(weighted_dd2, params.chart_sim_parameter,params.RHS, params.simulator.parallel, chart);
             [chart, ~, ~, ~, ~]                       = landmark(chart, params.RHS.t0, params.RHS.chi_p,params.RHS.threshold(1), params.RHS.connectivity_threshold);
              disp(['No. of landmarks after exploration is ', num2str(length(chart)) ])

             if params.relearn.option == 1
                  K_next                                   = length(chart);
                  index_to_learn                           = K+1: K_next;
                 [ chart ]                                 = relearn_chart(chart, params.relearn.settings, index_to_learn);
                 [chart, connectivity, P, bin_N, bins]     = landmark(chart, params.RHS.t0, params.RHS.chi_p,params.RHS.threshold(1), params.RHS.connectivity_threshold);
             end
             K                                         = length(chart);
             random_start;
             explore_round                             = explore_round + 1;
        end

        save( params.paths.chart ,'chart','connectivity','P')

        disp('exploration learning stage is completed')
