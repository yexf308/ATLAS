    Mean_FPT                           = zeros(2*3,10);
    relative_error_FPT                 = zeros(4,10);
    t_final_ori                        = zeros(2,1);
    t_final                            = zeros(4,10);

    
    theta_s                           = atan(b2^2/(2*b1));
    well_threshold                    = [-2*atan(4-sqrt(15)), -2*atan(4+sqrt(15))+2*pi];
    [FPT_ori1, t_final_ori(1)]        = MFPT_halfmoon_ori(well_threshold, output, phi_ori,N_IC,RHS_parameter,chart_sim_parameter, nonlin_trans_inv);

    well_threshold                    = [-2*atan(4+sqrt(15)), -2*atan(4-sqrt(15))]; 
    [FPT_ori2, t_final_ori(2)]        = MFPT_halfmoon_ori(well_threshold, output, phi_ori,N_IC,RHS_parameter,chart_sim_parameter, nonlin_trans_inv);


    %% ATLAS simulator
    for k = 1:10
        disp([num2str(k),'th round'])
        chart_fileName                     = [datapath,'chart',num2str(k),'.mat'];
        FPT_fileName                       = [datapath,'Halfmoon_FPT',num2str(k),'.mat'];
        load( chart_fileName )

        K                                  = length(chart);
        phi_store                          = zeros(1,K);
        for i = 1:K
             X_int                         = chart{i}.X_int;
             phi_store(i)                  = atan2(X_int(2), X_int(1));   
        end


        [chart, connectivity, P]           = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
        chart_sim_parameter.connectivity   = connectivity;  
        chart_sim_parameter.X_int          = chart{chart_sim_parameter.nearest}.X_int;
        mode                               = 2;


        well_threshold                     = [-2*atan(4-sqrt(15)), -2*atan(4+sqrt(15))+2*pi];
        disp('Fast mode projection')
        option                             = 1;
        weighted_dd2                       = @(X0, chart,neigh,nearest, connectivity_indices) weighted_drift_diffusion2( X0, chart, neigh, nearest, connectivity_indices, t0, chi_p, D,d, threshold, option,mode );
        [FPT_sim1, t_final(1,k)]                         = MFPT_halfmoon_sim(weighted_dd2, well_threshold, output, phi_ori, N_IC, phi_store, chart,chart_sim_parameter);

        disp('Orthogonal projection')
        option                             = 2;
        weighted_dd2                       = @(X0, chart,neigh,nearest, connectivity_indices) weighted_drift_diffusion2( X0, chart, neigh, nearest, connectivity_indices, t0, chi_p, D,d, threshold, option,mode );
        [FPT_ort1, t_final(2,k)]                         = MFPT_halfmoon_sim(weighted_dd2, well_threshold, output, phi_ori, N_IC, phi_store, chart,chart_sim_parameter);

        well_threshold                     =[-2*atan(4+sqrt(15)), -2*atan(4-sqrt(15))];
        disp('Fast mode projection')
        option                             = 1;
        weighted_dd2                       = @(X0, chart,neigh,nearest, connectivity_indices) weighted_drift_diffusion2( X0, chart, neigh, nearest, connectivity_indices, t0, chi_p, D,d, threshold, option,mode );
        [FPT_sim2, t_final(3,k)]                         = MFPT_halfmoon_sim(weighted_dd2, well_threshold, output, phi_ori, N_IC, phi_store, chart,chart_sim_parameter);

        disp('Orthogonal projection')
        option                             = 2;
        weighted_dd2                       = @(X0, chart,neigh,nearest, connectivity_indices) weighted_drift_diffusion2( X0, chart, neigh, nearest, connectivity_indices, t0, chi_p, D,d, threshold, option,mode );
       [FPT_ort2, t_final(4,k)]                          = MFPT_halfmoon_sim(weighted_dd2, well_threshold, output, phi_ori, N_IC, phi_store, chart, chart_sim_parameter);

       save(FPT_fileName,'FPT_ori1','FPT_ori2','FPT_sim1','FPT_sim2','FPT_ort1','FPT_ort2')

        Mean_FPT(1,k) = mean(FPT_ori1);
        Mean_FPT(2,k) = mean(FPT_ori2);
        Mean_FPT(3,k) = mean(FPT_sim1);
        Mean_FPT(4,k) = mean(FPT_sim2);
        Mean_FPT(5,k) = mean(FPT_ort1);
        Mean_FPT(6,k) = mean(FPT_ort2);

       if k ==1
          disp(['Original simulator'])
          disp(['The original simulator at right well has ',num2str(mean(FPT_ori1),'%5.0f'),'+-',num2str(std(FPT_ori1)*1.96/sqrt(N_IC),'%5.0f')])
          disp(['The original simulator at left well has ', num2str(mean(FPT_ori2),'%5.0f'),'+-',num2str(std(FPT_ori2)*1.96/sqrt(N_IC),'%5.0f')])
          disp(['Oblique projection'])
          disp(['The ATLAS simulator at right well has ',num2str(mean(FPT_sim1),'%5.0f'),'+-',num2str(std(FPT_sim1)*1.96/sqrt(N_IC),'%5.0f')])
          disp(['The ATLAS simulator at left well has ', num2str(mean(FPT_sim2),'%5.0f'),'+-',num2str(std(FPT_sim2)*1.96/sqrt(N_IC),'%5.0f')])
          disp(['Orthogonal projection'])
          disp(['The ATLAS simulator at right well has ',num2str(mean(FPT_ort1),'%5.0f'),'+-',num2str(std(FPT_ort1)*1.96/sqrt(N_IC),'%5.0f')])
          disp(['The ATLAS simulator at left well has ', num2str(mean(FPT_ort2),'%5.0f'),'+-',num2str(std(FPT_ort2)*1.96/sqrt(N_IC),'%5.0f')])

       end


        relative_error_FPT(1, k) = (Mean_FPT(3,k)-Mean_FPT(1,k))/Mean_FPT(1,k);
        relative_error_FPT(2, k) = (Mean_FPT(4,k)-Mean_FPT(2,k))/Mean_FPT(2,k);
        relative_error_FPT(3, k) = (Mean_FPT(5,k)-Mean_FPT(1,k))/Mean_FPT(1,k);
        relative_error_FPT(4, k) = (Mean_FPT(6,k)-Mean_FPT(2,k))/Mean_FPT(2,k);



    end

    disp(['Average relative error of Oblique ATLAS simulator at right state is ',    num2str(mean(relative_error_FPT(1,:)),'%5.4f'),'+-',num2str(std(relative_error_FPT(1,:))*1.96/sqrt(10),'%5.4f')])
    disp(['Average relative error at Oblique ATLAS simulator at left state  is ',    num2str(mean(relative_error_FPT(2,:)),'%5.4f'),'+-',num2str(std(relative_error_FPT(2,:))*1.96/sqrt(10),'%5.4f')])
    disp(['Average relative error of Orthogonal ATLAS simulator at right state is ', num2str(mean(relative_error_FPT(3,:)),'%5.4f'),'+-',num2str(std(relative_error_FPT(3,:))*1.96/sqrt(10),'%5.4f')])
    disp(['Average relative error at Orthogonal ATLAS simulator at left state is ',  num2str(mean(relative_error_FPT(4,:)),'%5.4f'),'+-',num2str(std(relative_error_FPT(4,:))*1.96/sqrt(10),'%5.4f')])


    save(Halfmoon_analysis_fileName,'Mean_FPT','relative_error_FPT','t_final_ori','t_final')