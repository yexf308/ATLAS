
       Mean_FPT                           = zeros(6,10);
       relative_error_FPT                 = zeros(3,10);
       Q1                                 = zeros(3,3,10);
       Q2                                 = zeros(3,3,10);
       t_final_ori                        = zeros(3,1);
       t_final                            = zeros(3,10);
       dt_s                               = t0;
    
        disp('Starting MFPT part')
        well_threshold1       = [-pi/3, pi/3]; 
        well_threshold2       = [pi/3, pi];
        well_threshold3       = [-pi, -pi/3];


        disp('Original simulator')
        [FPT_ori1,exit_curr1_ori, t_final_ori(1)] = MFPT_butane2_ori(well_threshold1, output, N_IC, RHS_parameter,chart_sim_parameter);
        [FPT_ori2,exit_curr2_ori, t_final_ori(2)] = MFPT_butane2_ori(well_threshold2, output, N_IC, RHS_parameter,chart_sim_parameter);
        [FPT_ori3,exit_curr3_ori, t_final_ori(3)] = MFPT_butane2_ori(well_threshold3, output, N_IC, RHS_parameter,chart_sim_parameter);

    for k = 1:10
        disp([num2str(k),'th round'])
        chart_fileName                     = [datapath,'chart',num2str(k),'.mat'];
        FPT_fileName                       = [datapath,'Butane_FPT',num2str(k),'.mat'];
        load( chart_fileName )

        K                                  = length(chart);
        phi_store                          = zeros(1,K);
        for i = 1:K
             X_int                         = chart{i}.X_int;
             phi_store(i)                  = atan2(X_int(6), X_int(4));   
        end

        [chart, connectivity, P]           = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
        chart_sim_parameter.connectivity   = connectivity;  
        chart_sim_parameter.X_int          = chart{chart_sim_parameter.nearest}.X_int;
        mode                               = 2;

        option                             = 1;
        weighted_dd2                       = @(X0, chart,neigh,nearest, connectivity_indices) weighted_drift_diffusion2( X0, chart, neigh, nearest, connectivity_indices, t0, chi_p, D,d, threshold, option,mode );

        [FPT_sim1,exit_curr1_sim,t_final(1,k)] = MFPT_butane2_sim(weighted_dd2, well_threshold1, output, N_IC, phi_store, chart,chart_sim_parameter);
        [FPT_sim2,exit_curr2_sim,t_final(2,k)] = MFPT_butane2_sim(weighted_dd2, well_threshold2, output, N_IC, phi_store, chart,chart_sim_parameter);
        [FPT_sim3,exit_curr3_sim,t_final(1,k)] = MFPT_butane2_sim(weighted_dd2, well_threshold3, output, N_IC, phi_store, chart,chart_sim_parameter);


        save(FPT_fileName, 'FPT_ori1', 'FPT_ori2', 'FPT_ori3','exit_curr1_ori', 'exit_curr2_ori', 'exit_curr3_ori', ...
                           'FPT_sim1', 'FPT_sim2', 'FPT_sim3','exit_curr1_sim', 'exit_curr2_sim', 'exit_curr3_sim')

        Mean_FPT(1,k) = mean(FPT_ori1);
        Mean_FPT(2,k) = mean(FPT_ori2);
        Mean_FPT(3,k) = mean(FPT_ori3);
        Mean_FPT(4,k) = mean(FPT_sim1);
        Mean_FPT(5,k) = mean(FPT_sim2);
        Mean_FPT(6,k) = mean(FPT_sim3);

        relative_error_FPT(1, k) = (Mean_FPT(4,k)-Mean_FPT(1,k))/Mean_FPT(1,k);
        relative_error_FPT(2, k) = (Mean_FPT(5,k)-Mean_FPT(2,k))/Mean_FPT(2,k);
        relative_error_FPT(3, k) = (Mean_FPT(6,k)-Mean_FPT(3,k))/Mean_FPT(3,k);   


         Q1(1,1,k) = -1/mean(FPT_ori1);
         Q1(2,2,k) = -1/mean(FPT_ori2);
         Q1(3,3,k) = -1/mean(FPT_ori3);
         Q1(1,2,k) = sum(exit_curr1_ori>0)/N_IC * (-Q1(1,1,k));
         Q1(1,3,k) = sum(exit_curr1_ori<0)/N_IC * (-Q1(1,1,k));
         Q1(2,1,k) = sum(exit_curr2_ori>0)/N_IC * (-Q1(2,2,k));
         Q1(2,3,k) = sum(exit_curr2_ori<0)/N_IC * (-Q1(2,2,k));
         Q1(3,1,k) = sum(exit_curr3_ori<0)/N_IC * (-Q1(3,3,k));
         Q1(3,2,k) = sum(exit_curr3_ori>0)/N_IC * (-Q1(3,3,k));


         Q2(1,1,k) = -1/mean(FPT_sim1);
         Q2(2,2,k) = -1/mean(FPT_sim2);
         Q2(3,3,k) = -1/mean(FPT_sim3);
         Q2(1,2,k) = sum(exit_curr1_sim>0)/N_IC * (-Q2(1,1,k));
         Q2(1,3,k) = sum(exit_curr1_sim<0)/N_IC * (-Q2(1,1,k));
         Q2(2,1,k) = sum(exit_curr2_sim>0)/N_IC * (-Q2(2,2,k));
         Q2(2,3,k) = sum(exit_curr2_sim<0)/N_IC * (-Q2(2,2,k));
         Q2(3,1,k) = sum(exit_curr3_sim<0)/N_IC * (-Q2(3,3,k));
         Q2(3,2,k) = sum(exit_curr3_sim>0)/N_IC * (-Q2(3,3,k));


        if k ==1
            disp(['The original simulator at trans state has ',      num2str(mean(FPT_ori1),'%5.4f'),'+-',num2str(std(FPT_ori1)*1.96/sqrt(N_IC),'%5.4f')])
            disp(['The original simulator at cis-top state has ',    num2str(mean(FPT_ori2),'%5.4f'),'+-',num2str(std(FPT_ori2)*1.96/sqrt(N_IC),'%5.4f')])
            disp(['The original simulator at cis-bottom state has ', num2str(mean(FPT_ori3),'%5.4f'),'+-',num2str(std(FPT_ori3)*1.96/sqrt(N_IC),'%5.4f')])
            disp(['The ATLAS simulator at trans state has ',         num2str(mean(FPT_sim1),'%5.4f'),'+-',num2str(std(FPT_sim1)*1.96/sqrt(N_IC),'%5.4f')])
            disp(['The ATLAS simulator at cis-top state has ',       num2str(mean(FPT_sim2),'%5.4f'),'+-',num2str(std(FPT_sim2)*1.96/sqrt(N_IC),'%5.4f')])
            disp(['The ATLAS simulator at cis-bottom state has ',    num2str(mean(FPT_sim3),'%5.4f'),'+-',num2str(std(FPT_sim3)*1.96/sqrt(N_IC),'%5.4f')])

        end


    end
    disp(['Average relative error at Trans state is ',      num2str(mean(relative_error_FPT(1,:)),'%5.4f'),'+-',num2str(std(relative_error_FPT(1,:))*1.96/sqrt(10),'%5.4f')])
    disp(['Average relative error at cis-top state is ',    num2str(mean(relative_error_FPT(2,:)),'%5.4f'),'+-',num2str(std(relative_error_FPT(2,:))*1.96/sqrt(10),'%5.4f')])
    disp(['Average relative error at cis-bottom state is ', num2str(mean(relative_error_FPT(3,:)),'%5.4f'),'+-',num2str(std(relative_error_FPT(3,:))*1.96/sqrt(10),'%5.4f')])

    disp(['The mean Q1(1,1) ', num2str(mean(Q1(1,1,:)),'%5.4f'),'+-', num2str(std(Q1(1,1,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q1(1,2) ', num2str(mean(Q1(1,2,:)),'%5.4f'),'+-', num2str(std(Q1(1,2,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q1(1,3) ', num2str(mean(Q1(1,3,:)),'%5.4f'),'+-', num2str(std(Q1(1,3,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q1(2,1) ', num2str(mean(Q1(2,1,:)),'%5.4f'),'+-', num2str(std(Q1(2,1,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q1(2,2) ', num2str(mean(Q1(2,2,:)),'%5.4f'),'+-', num2str(std(Q1(2,2,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q1(2,3) ', num2str(mean(Q1(2,3,:)),'%5.4f'),'+-', num2str(std(Q1(2,3,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q1(3,1) ', num2str(mean(Q1(3,1,:)),'%5.4f'),'+-', num2str(std(Q1(3,1,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q1(3,2) ', num2str(mean(Q1(3,2,:)),'%5.4f'),'+-', num2str(std(Q1(3,2,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q1(3,3) ', num2str(mean(Q1(3,3,:)),'%5.4f'),'+-', num2str(std(Q1(3,3,:))*1.96/sqrt(10),'%5.4f')])

    disp(['The mean Q2(1,1) ', num2str(mean(Q2(1,1,:)),'%5.4f'),'+-', num2str(std(Q2(1,1,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q2(1,2) ', num2str(mean(Q2(1,2,:)),'%5.4f'),'+-', num2str(std(Q2(1,2,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q2(1,3) ', num2str(mean(Q2(1,3,:)),'%5.4f'),'+-', num2str(std(Q2(1,3,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q2(2,1) ', num2str(mean(Q2(2,1,:)),'%5.4f'),'+-', num2str(std(Q2(2,1,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q2(2,2) ', num2str(mean(Q2(2,2,:)),'%5.4f'),'+-', num2str(std(Q2(2,2,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q2(2,3) ', num2str(mean(Q2(2,3,:)),'%5.4f'),'+-', num2str(std(Q2(2,3,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q2(3,1) ', num2str(mean(Q2(3,1,:)),'%5.4f'),'+-', num2str(std(Q2(3,1,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q2(3,2) ', num2str(mean(Q2(3,2,:)),'%5.4f'),'+-', num2str(std(Q2(3,2,:))*1.96/sqrt(10),'%5.4f')])
    disp(['The mean Q2(3,3) ', num2str(mean(Q2(3,3,:)),'%5.4f'),'+-', num2str(std(Q2(3,3,:))*1.96/sqrt(10),'%5.4f')])
    save(Butane_analysis_fileName,'Mean_FPT','relative_error_FPT','Q1','Q2','t_final')
    