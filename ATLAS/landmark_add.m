function   [ chart, connectivity,  P]   =  landmark_add(chart, connectivity, chart_add, t0, chi_p,threshold, connectivity_threshold)
K                          = length(chart);
connectivity_new           = zeros(K+1,K+1);
connectivity_new(1:K, 1:K) = connectivity;
chart                      = [chart chart_add];

for i=1:K+1
     if    check_x_in_chart( chart{i}.X_int, chart{K+1}, t0, connectivity_threshold, chi_p,threshold ) ...
           || check_x_in_chart( chart{K+1}.X_int, chart{i}, t0, connectivity_threshold, chi_p,threshold ) 
              dist1   =  distance_chart_to_x(chart{i}.X_int, chart{K+1}, t0, chi_p,threshold);
              dist2   =  distance_chart_to_x(chart{K+1}.X_int, chart{i}, t0, chi_p,threshold);
                  connectivity_new(K+1,i) = max(dist1,dist2);
       
     end
end
connectivity_new(K+1, K+1) = 1;
connectivity_new(:, K+1)   = connectivity_new(K+1,:)';
connectivity      = connectivity_new;
P                 = graph(connectivity,'upper');
end

