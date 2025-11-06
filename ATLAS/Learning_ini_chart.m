% Find the chart for randomly sampled initial points on the manifold.
% These initial points are evenly sampled from a long trajectory.

chart                   = cell(1, params.learning.K_int);
Sim_parameter           = struct(                        ...
                                'T_max',        params.learning.T_one,   ...
                                'dt',           params.RHS.dt ,          ...
                                'gap',          1,                       ...
                                'X_int',        params.learning.X_int    ...
                         );


disp(['Initial point: ', num2str(params.learning.X_int), ' T_max: ', num2str(Sim_parameter.T_max)])

% Simulate a long trajectory starting from this initial condition.
tic

[output]               = params.simulator.one(Sim_parameter);
t1=toc;
disp(['one long trajectory is simulated. The time spent is ', num2str(t1/60), ' mins.' ])



% evenly sample N initial points from the trajectory
output                 = output(10000:end,:);
[L, ~]                 = size(output);
X_int_sample           = output(1: round(L/params.learning.K_int):end, :);

tic
disp(['Parallelly run chart simulator starting from ', num2str(params.learning.K_int) ,' initial points. '])
for i = 1 : params.learning.K_int
       disp(['Landmark',num2str(i), ' is learned.'])

       Sim_parameter                 = struct(                           ...
                                           'T_max',        params.RHS.T_max,        ...
                                           'dt',           params.RHS.dt ,          ...
                                           'N',            params.RHS.N,            ...
                                           'X_int',        X_int_sample(i,:)     ...
                                      );
       [data, Cov_store, Mean_store] = params.simulator.parallel(Sim_parameter);
       [chart{i}]                    = Learning_Slow_Manifold( data, Cov_store, Mean_store, params.RHS.D, params.RHS.d, params.RHS.modify );


end
t2 =toc;
disp(['Initial learning stage is completed. The time spent is ', num2str(t2/60), ' mins.' ])
  
