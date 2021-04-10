% previous relearning is to make sure it relaxes onto the manifold. 
% Now we set N 10 times bigger.  
% The upperbound and lowerbound are set to [10dt, 15dt]. This is the
% relaxation timescale. 


         RHS_parameter_relearn.T_max      = 15*dt;
         RHS_parameter_relearn.UpperBound = 15*dt;
         RHS_parameter_relearn.LowerBound = 10*dt;
         RHS_parameter_relearn.modify     = 1;
         relearn_parameter      = struct(                            ...
                    'N_relearn',          10*N,                      ...
                    'iter',               2,                         ...
                    'RHS_parameter',      RHS_parameter_relearn,     ...
                    'relative_threshold',  [0.005, 0.02],            ...
                    'simulator_par',      @(Sim_parameter) simEuler_par(RHS_parameter_relearn,Sim_parameter) ...
                    );       

        disp('Starting  Relearning stage')
        K                                     = length(chart);        
        index_to_learn                        = 1:K;
        [ chart,~ ]                           = relearn_chart(chart, relearn_parameter, index_to_learn);
        [chart, connectivity, P]              = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
        save( chart_fileName ,'chart','connectivity','P')
               
        disp('Relearning stage is completed.')
     