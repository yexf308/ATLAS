
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
        chart_fileName                     = [datapath,'chart',num2str(k),'.mat'];
        chart_part_fileName                = [datapath,'chart_part',num2str(k),'.mat'];
        TranM_fileName                     = [datapath,'TranM',num2str(k),'.mat'];
        diary_fileName                     = [datapath,'diary',num2str(k),'.txt'];
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
    diary_FPT_fileName                     = [datapath,'diary_FPT.txt'];
    Butane_analysis_fileName               = [datapath,'Butane_analysis.mat'];
    diary(diary_FPT_fileName)

     
    %% Simulate long trajectory with the original simulator
    
    butane_ori_simulation;

    %% Calculate MFPT
    
    butane_MFPT;
    
    diary off

end