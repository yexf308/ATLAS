      %% Generate the chart for randomly sampled initial points.

        disp('Starting initial learning stage')
        Learning_ini_chart;
        [chart, ~, ~, ~]                        = landmark(chart, params.RHS.t0, params.RHS.chi_p,params.RHS.threshold(1), params.RHS.connectivity_threshold);

        if params.relearn.option == 1
            K                                       = length(chart);
            index_to_learn                          = 1:K;
            [ chart ]                               = relearn_chart(chart, params.relearn.settings, index_to_learn);
            [chart, connectivity, P, bin_N, bins]   = landmark(chart, params.RHS.t0, params.RHS.chi_p,params.RHS.threshold(1), params.RHS.connectivity_threshold);
        end
        save(params.paths.chart_part , 'chart','connectivity','P');
        disp(['No. of landmarks after initial learning stage is ', num2str(length(chart)) ])

        disp('Initial learning stage is completed')
