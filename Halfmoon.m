
clear all
close all
clc


system_type        = 'Halfmoon';
disp(['Starting ',system_type,' model.'])
addpath_settings;

prompt = '1-Learning ATLAS, 2-Simulating MRT';
x = input(prompt);

if x==1

    for k =1:10
    %% Set up the parameter
        rng(k)    
        set_parameter;
        disp('Parameter is set')  

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

         %% Using ATLAS simulator to run lots of trajs to further learn.

        halfmoon_simulation_learning;
           
        %% Relearn ATLAS    
        
        halfmoon_relearning;

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
    Halfmoon_analysis_fileName             = [datapath,'Halfmoon_analysis.mat'];
    diary(diary_FPT_fileName)


   %% Simulate long trajectory with the original simulator
   
    halfmoon_ori_simulation;
   
   %% Calculate MFPT 
    
    halfmoon_MFPT;
    
    diary off 


end