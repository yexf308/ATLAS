        RHS_parameter_relearn.modify     = 1;
        relearn_parameter      = struct(                            ...
                    'N_relearn',          2*N,                       ...
                    'iter',               1,                         ...
                    'RHS_parameter',      RHS_parameter_relearn,     ...
                    'relative_threshold',  [0.005, 0.02],            ...
                    'simulator_par',      @(Sim_parameter) simEulernonlintrans_par(RHS_parameter_relearn,Sim_parameter) ...
                    );   


        disp('Starting Relearning Chart')
        K                                     = length(chart);        
        index_to_learn                        = 1:K;
        [ chart,~ ]                           = relearn_chart(chart, relearn_parameter, index_to_learn);
        [chart, connectivity, P]              = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
        save( chart_fileName ,'chart','connectivity','P')  
        
        disp('Relearning stage is completed.')
