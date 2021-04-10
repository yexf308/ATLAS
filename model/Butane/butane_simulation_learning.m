       disp('Use ATLAS simulator to run long trajs')

        for i = 1:500
            disp(['traj No.',num2str(i) ])
            K                                        = length(chart);
            chart_sim_parameter.nearest              = max(round(rand * K),1);
            chart_sim_parameter.connectivity         = connectivity;  
            chart_sim_parameter.X_int                = chart{chart_sim_parameter.nearest}.X_int;
            chart_sim_parameter.Nstep                = 2*10^3;
            [~, ~, chart]                            = ATLAS_simulator2(weighted_dd2, chart_sim_parameter,RHS_parameter, simulator_par, chart);        
            [chart, connectivity, ~]                 = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
            if relearn_option == 1
                K_next                                   = length(chart); 
                if K_next>K
                     index_to_learn                            = K+1: K_next;
                     [ chart ]                                 = relearn_chart(chart, relearn_parameter, index_to_learn);
                     [chart, connectivity, P, bin_N ]          = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
                     save( chart_fileName ,'chart','connectivity','P')
                end
            end
            
            save( chart_fileName ,'chart','connectivity','P')       
        end 
        disp('Simulation stage is completed')
