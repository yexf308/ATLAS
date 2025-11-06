
   disp('Starting MSM learning stage')

   [TranM]                                           = MSM(chart, connectivity, params.MSM, weighted_dd2, params.RHS.d);
    step                                             = params.MSM.step;
    N_state                                          = params.MSM.N_state;
    dt_s                                             = params.MSM.dt_s;
    save( params.paths.TranM ,'TranM','step','N_state','dt_s')
    disp('MSM learning stage is completed')