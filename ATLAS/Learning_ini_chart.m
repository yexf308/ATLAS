% Find the chart for randomly sampled initial points on the manifold. 
%These initial points are evenly sampled from a long trajectory
 
chart                   = cell(1,K_int);   
Sim_parameter           = struct(                        ...
                                'T_max',        T_one,   ...
                                'dt',           dt ,     ...
                                'gap',          1,       ...
                                'X_int',        X_int    ...
                         );
                 
                     
disp(['Initial point: ', num2str(X_int), ' T_max: ', num2str(Sim_parameter.T_max)])   

% Simulate a long trajctory starting from this initial condition. 
tic 
  
[output]               = simulator_one(Sim_parameter);  
t1=toc;  
disp(['one long trajctory is simulated. The time spent is ', num2str(t1/60), ' mins.' ])
 
 

% evenly sample N initial points from the trajctory
output                 = output(10000:end,:);
[L, ~]                 = size(output); 
X_int_sample           = output(1: round(L/K_int):end, :);

tic
disp(['Parallelly run chart simulator starting from ', num2str(K_int) ,' initial points. '])
for i = 1 : K_int
       disp(['Landmark',num2str(i), ' is learned.'])

       Sim_parameter                 = struct(                           ...
                                           'T_max',        T_max,        ...
                                           'dt',           dt ,          ...
                                           'N',            N,            ...
                                           'X_int',        X_int_sample(i,:)     ...
                                      );
       [data, Cov_store, Mean_store] = simulator_par(Sim_parameter); 
       [chart{i}]                    = Learning_Slow_Manifold( data, Cov_store, Mean_store, D, d,modify );
      
      
end
t2 =toc;
disp(['Initial learning stage is completed. The time spent is ', num2str(t2/60), ' mins.' ])
  
