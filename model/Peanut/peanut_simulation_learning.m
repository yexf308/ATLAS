      disp('Use ATLAS simulator to simulate a single traj.')
      chart_sim_parameter.Nstep                 = 2*10^7;
      chart_sim_parameter.gap                   = 1000; 
      [~, ~, chart]                             = ATLAS_simulator2(weighted_dd2, chart_sim_parameter, RHS_parameter, simulator_par, chart);
      [chart, connectivity, P]                  = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);  
      save( chart_fileName ,'chart','connectivity','P')