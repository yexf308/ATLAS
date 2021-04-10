      %% Generate the chart for randomly sampled initial points.

        disp('Starting initial learning stage')
        Learning_ini_chart;
        [chart, ~, ~, ~]                        = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
        
        if relearn_option == 1
            K                                       = length(chart);
            index_to_learn                          = 1:K;
            [ chart ]                               = relearn_chart(chart, relearn_parameter, index_to_learn);
            [chart, connectivity, P, bin_N, bins]   = landmark(chart, t0, chi_p,threshold(1), connectivity_threshold);
        end
        save(chart_part_fileName , 'chart','connectivity','P'); 
        disp(['No. of landmarks after initial learning stage is ', num2str(length(chart)) ])

        disp('Initial learning stage is completed')
