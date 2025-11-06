
clear all 
close all
clc


system_type        = 'Butane'; 
disp(['Starting ',system_type,' model.'])
addpath_settings; % setup path

prompt = '1-Learning ATLAS, 2-Simulating MRT';
x = input(prompt);

if x==1

    for k =1:10

        %% Set up the parameter
        rng(k) 
        set_parameter;    % setup parameter
        disp('Parameter is set') 

        %set path 
        butane_paths                       = build_butane_paths(datapath, k);
        chart_fileName                     = butane_paths.chart;
        chart_part_fileName                = butane_paths.chart_part;
        TranM_fileName                     = butane_paths.TranM;
        diary_fileName                     = butane_paths.diary;
        diary(diary_fileName)

        tstart=tic; 

        %% Generate the chart for randomly sampled initial points.

        initial_learning;


        %% Use the exploration mode to further Learn 
         
        exploration_learning;

        %% Using ATLAS simulator to run lots of trajs to further Learn

        butane_simulation_learning;


        %% Relearn ATLAS 
        
        butane_relearning;

        %% Build up Markov state model

        MSM_learning;
        
        
        %% 
        t_final      = toc(tstart);
        disp(['Learning stage is finished. The time spent is ',num2str(t_final/3600), ' hours']);
        diary off

    end
    
elseif x==2

    set_parameter;
    disp('Parameter is set') 
    butane_paths                           = build_butane_paths(datapath);
    diary_FPT_fileName                     = butane_paths.diary_FPT;
    Butane_analysis_fileName               = butane_paths.Butane_analysis;
    diary(diary_FPT_fileName)

     
    %% Simulate long trajectory with the original simulator
    
    butane_ori_simulation;

    %% Calculate MFPT
    
    butane_MFPT;
    
    diary off

end