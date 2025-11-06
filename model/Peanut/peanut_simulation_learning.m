      disp('Use ATLAS simulator to simulate a single traj.')
      params.chart_sim_parameter.Nstep                 = 2*10^7;
      params.chart_sim_parameter.gap                   = 1000;
      [~, ~, chart]                             = ATLAS_simulator2(weighted_dd2, params.chart_sim_parameter, params.RHS, params.simulator.parallel, chart);
      [chart, connectivity, P]                  = landmark(chart, params.RHS.t0, params.RHS.chi_p,params.RHS.threshold(1), params.RHS.connectivity_threshold);
      save( params.paths.chart ,'chart','connectivity','P')
