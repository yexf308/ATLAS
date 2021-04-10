%% Relearn the chart with different diffusion coefficient or to further relax.
% The diffusion coefficient here is smaller and exploration mode is not
% needed since the invariant manifold is the same. 

% another situation is when T_max is not enough to relax the IC onto the
% invariant manifold because it might be too far away or different
% compoents have different relaxation scale.




%set up new parameter if needed.
%note t0, lowerbound and upperbound of training window may be reset. 

function [ chart, i_to_learn ] = relearn_chart(chart, relearn_parameter, index_to_learn)

N_relearn          = relearn_parameter.N_relearn;
RHS_parameter      = relearn_parameter.RHS_parameter;
simulator_par      = relearn_parameter.simulator_par;
iter_max           = relearn_parameter.iter;
relative_threshold = relearn_parameter.relative_threshold;
modify             = RHS_parameter.modify;
T_max              = RHS_parameter.T_max;
dt                 = RHS_parameter.dt;
D                  = RHS_parameter.D;
d                  = RHS_parameter.d;
t0                 = RHS_parameter.t0;
chi_p              = RHS_parameter.chi_p;
threshold          = RHS_parameter.threshold;
connectivity_threshold = RHS_parameter.connectivity_threshold;

K             = length(chart);

indicator = ones(1,K);  % indicate which landmark should keep 
i_to_learn = [];
for i = 1:K
     
     if ismember(i, index_to_learn)
         disp(['Landmark No. ',num2str(i)])
         relative_error  =  ones(1,2);
         j               = 1;  % No. of iter 
         while (relative_error(1)>relative_threshold(1) || relative_error(2)>relative_threshold(2) ) && indicator(i) == 1 && j<=iter_max
             % if both relative error are less than 2%, or landmark is close
             % to other previous landmarks or reaches maximum, jump out of
             % loop
             
            
            disp(['iter No. ',num2str(j)])
            Sim_parameter                    = struct(                               ...
                                               'T_max',        T_max,                ...
                                               'dt',           dt ,                  ...
                                               'N',            N_relearn,            ...
                                               'X_int',        chart{i}.X_int        ...
                                             );
            [data, Cov_store, Mean_store]       = simulator_par(Sim_parameter); 
            [chart_new]                         = Learning_Slow_Manifold( data, Cov_store, Mean_store, D, d,modify );
            relative_error(1)                   = norm(chart{i}.sigma - chart_new.sigma) / norm(chart_new.sigma);
            relative_error(2)                   = norm(chart{i}.b - chart_new.b)/norm(chart_new.b);
            disp(['relative error is ', num2str(relative_error)])
            chart{i}                            = chart_new; 
            
            
            for k = 1:i-1
                 if  indicator(k) == 1 && check_x_in_chart( chart{i}.X_int, chart{k}, t0, 0.3, chi_p,threshold(1) ) && check_x_in_chart( chart{k}.X_int, chart{i}, t0, 0.3, chi_p,threshold(1) )
                     indicator(i) = 0 ; % this chart is been removed.
                     disp(['Landmark ', num2str(i), 'is removed'])
                     break
                 end
            end
            
            
            j = j+1;     
       
            
         end
         
          if j>iter_max && indicator(i) == 1
              i_to_learn = [i_to_learn, i]; % these landmarks are not accuate 
          end
         
     end
     
    
end

chart                 = chart(find(indicator));
[chart, ~, ~]         = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);


