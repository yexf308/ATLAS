% previous relearning is to make sure it relaxes onto the manifold.
% Now we set N 10 times bigger.
% The upper bound and lower bound are set to [10dt, 15dt]. This is the
% relaxation timescale.


         params.relearn.RHS.T_max      = 15*params.RHS.dt;
         params.relearn.RHS.UpperBound = 15*params.RHS.dt;
         params.relearn.RHS.LowerBound = 10*params.RHS.dt;
         params.relearn.RHS.modify     = 1;
         relearn_parameter      = struct(                            ...
                    'N_relearn',          10*params.RHS.N,          ...
                    'iter',               2,                         ...
                    'RHS_parameter',      params.relearn.RHS,       ...
                    'relative_threshold',  [0.005, 0.02],            ...
                    'simulator_par',      @(Sim_parameter) simEuler_par(params.relearn.RHS,Sim_parameter) ...
                    );

        disp('Starting  Relearning stage')
        K                                     = length(chart);
        index_to_learn                        = 1:K;
        [ chart,~ ]                           = relearn_chart(chart, relearn_parameter, index_to_learn);
        [chart, connectivity, P]              = landmark(chart, params.RHS.t0, params.RHS.chi_p,params.RHS.threshold(1), params.RHS.connectivity_threshold);
        save( params.paths.chart ,'chart','connectivity','P')

        disp('Relearning stage is completed.')
     
