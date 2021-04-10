
clear all
close all
clc

for k =1:10
    
rng(k) 

%% Set up the parameter

system_type        = 'Peanut'; 
relearn_option     = 0;
disp(['Starting ',system_type,' model.'])

addpath_settings;
set_parameter;
disp('Parameter is set') 


chart_fileName                     = [datapath,'chart',num2str(k),'.mat'];
chart_part_fileName                = [datapath,'chart_part',num2str(k),'.mat'];
TranM_fileName                     = [datapath,'TranM',num2str(k),'.mat'];
FPT_fileName                       = [datapath,'Peanut_FPT',num2str(k),'.mat'];
diary_fileName                     = [datapath,'diary',num2str(k),'.txt'];
diary(diary_fileName)
tstart=tic; 
    
%% Generate the chart for randomly sampled initial points on the manifold.
    disp('Starting initial learning stage')
    Learning_ini_chart;
    [chart, connectivity, P, bin_N, bins]        = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
    disp(['No. of landmarks is ', num2str(length(chart)) ])
    save(chart_part_fileName , 'chart','connectivity','P');  

    %% weighted_dd set to mode 1: exploration mode

    mode                    = 1;
    weighted_dd2            = @(X0, chart,neigh,nearest, connectivity_indices) weighted_drift_diffusion2( X0, chart, neigh, nearest, connectivity_indices, t0, chi_p, D,d, threshold, option,mode );

    %% Use the exploration mode to further Learn 
    explore_round = 1;
    random_start;  
    while bin_N~=1 % if graph is not fully connected. This doesn't guarantee neighbour of landmarks cover the invariant manifold
          disp(['Round ',num2str(explore_round),' of exploration'])
          
          [X, nearest_store, chart]                      = ATLAS_simulator2(weighted_dd2, chart_sim_parameter,RHS_parameter, simulator_par, chart);      
          [chart, connectivity, P, bin_N, bins]          = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
          
          disp(['No. of landmarks after exploration is ', num2str(length(chart)) ])
          random_start;     
          explore_round                                  = explore_round + 1;
    end
    save( chart_fileName ,'chart','connectivity','P')
          
%%   Using ATLAS simulator to run a long traj.          
  chart_sim_parameter.Nstep                 = 2*10^7;
  chart_sim_parameter.gap                   = 1000; 
  [~, ~, chart]                             = ATLAS_simulator2(weighted_dd2, chart_sim_parameter, RHS_parameter, simulator_par, chart);
  [chart, connectivity, P]                  = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);  
  save( chart_fileName ,'chart','connectivity','P')


%% Build up Markov state model

    disp('Starting MSM part')
   [TranM]                                 = MSM(chart, connectivity, MSM_parameter, weighted_dd2, d);
    step                                   = MSM_parameter.step;
    N_state                                = MSM_parameter.N_state;
    dt_s                                   = MSM_parameter.dt_s;
    save( TranM_fileName ,'TranM','step','N_state','dt_s')  


    t_final      = toc(tstart);
    disp(['Learning stage is finished. The time spent is ',num2str(t_final/3600), ' hours']);


     diary off
end
 


