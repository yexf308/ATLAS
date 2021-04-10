% randomly pick an element in the smallest connected component in the
    % graph. Assign it to chart_sim_parameter. 
    bins_unique                              = unique(bins);
    bins_hist                                = histc(bins, bins_unique);
    [N_min, index_min]                       = min(bins_hist);
    minimum_block                            = find(bins==bins_unique(index_min));
    nearest                                  = minimum_block(ceil(rand * N_min));
    chart_sim_parameter.nearest              = nearest;
    chart_sim_parameter.connectivity         = connectivity;
    chart_sim_parameter.X_int                = chart{nearest}.X_int;
    
    
