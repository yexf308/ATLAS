
   disp('Starting MSM learning stage')

   [TranM]                                           = MSM(chart, connectivity, MSM_parameter, weighted_dd2, d);
    step                                             = MSM_parameter.step;
    N_state                                          = MSM_parameter.N_state;
    dt_s                                             = MSM_parameter.dt_s;
    save( TranM_fileName ,'TranM','step','N_state','dt_s')
    disp('MSM learning stage is completed')