 %% Original simulator
     disp('Use original simulation to simulate one single traj.')

    Sim_parameter                      = struct(        ...
                             'T_max',        T_one*3,   ...
                             'dt',           dt ,       ...
                             'gap',          1,         ...
                             'X_int',        X_int      ...
        );

    [output]                          = simulator_one(Sim_parameter);
    [L_out,~]                         = size(output);
    temp                              = nonlin_trans_inv(output');
    phi_ori                           = temp(1,:);
     disp(['one single traj is simulated with the original simulator'])
