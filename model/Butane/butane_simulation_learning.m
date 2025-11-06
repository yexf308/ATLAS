       disp('Use ATLAS simulator to run long trajs')

        for i = 1:500
            disp(['traj No.',num2str(i) ])
            K                                        = length(chart);
            params.chart_sim_parameter.nearest       = max(round(rand * K),1);
            params.chart_sim_parameter.connectivity  = connectivity;
            params.chart_sim_parameter.X_int         = chart{params.chart_sim_parameter.nearest}.X_int;
            params.chart_sim_parameter.Nstep         = 2*10^3;
            [~, ~, chart]                            = ATLAS_simulator2(weighted_dd2, params.chart_sim_parameter,params.RHS, params.simulator.parallel, chart);
            [chart, connectivity, ~]                 = landmark(chart, params.RHS.t0, params.RHS.chi_p,params.RHS.threshold(1), params.RHS.connectivity_threshold);
            if params.relearn.option == 1
                K_next                                   = length(chart);
                if K_next>K
                     index_to_learn                            = K+1: K_next;
                     [ chart ]                                 = relearn_chart(chart, params.relearn.settings, index_to_learn);
                     [chart, connectivity, P, bin_N ]          = landmark(chart, params.RHS.t0, params.RHS.chi_p,params.RHS.threshold(1), params.RHS.connectivity_threshold);
                     save( params.paths.chart ,'chart','connectivity','P')
                end
            end

            save( params.paths.chart ,'chart','connectivity','P')
        end
        disp('Simulation stage is completed')
