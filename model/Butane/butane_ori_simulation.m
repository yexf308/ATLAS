    disp('Use original simulation to simulate one single traj.')
    T_one                              = 100;
    N_IC                               = 30000;
    Sim_parameter           = struct(                        ...
        'T_max',        T_one*5,   ...
        'dt',           dt ,     ...
        'gap',          1,       ...
        'X_int',        X_int    ...
        );
    [output]                      = simulator_one(Sim_parameter);
    output                        = output(1000:end,:);
    phi_ori                       = atan2(output(:,6),output(:,4));
    disp(['one single traj is simulated with the original simulator'])
    