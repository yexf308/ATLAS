function [chart, connectivity, P, bin_N, bins] = landmark(chart_ori,t0, chi_p,threshold,connectivity_threshold)
% This function finds the landmark among the sample points and build the
% connectivity matrix for these landmark
% distance is diffusion distance with unit of sqrt(t0). 
K                = length(chart_ori);
indicator        = ones(1,K);


% if dist(i,j)<1-sqrt(1/2) and dist(j,i)<1-sqrt(1/2), remove j
for i   = 1:K
    for j= i+1:K
        if  indicator(i) &&  indicator(j)   ...
            && check_x_in_chart( chart_ori{i}.X_int, chart_ori{j}, t0, 0.3, chi_p,threshold ) ...
            && check_x_in_chart( chart_ori{j}.X_int, chart_ori{i}, t0, 0.3, chi_p ,threshold)   
        
                 indicator(j)       = 0;   
        end
    end
end

n                 = 1:K;
chart             = chart_ori( n(indicator>0) );



% i and j are connected if dist(i,j)<connectivity_threshold or dist(j,i)<connectivity_threshold
connectivity      = zeros(length(chart), length(chart));

for i   = 1:length(chart)
    connectivity(i,i) = 1;
    for j   = i+1:length(chart)
        if    check_x_in_chart( chart{i}.X_int, chart{j}, t0, connectivity_threshold, chi_p,threshold ) ...
           && check_x_in_chart( chart{j}.X_int, chart{i}, t0, connectivity_threshold, chi_p,threshold )  
                  dist1   =  distance_chart_to_x(chart{i}.X_int, chart{j}, t0, chi_p,threshold);
                  dist2   =  distance_chart_to_x(chart{j}.X_int, chart{i}, t0, chi_p,threshold);
                
                  connectivity(i,j) = max(dist1,dist2);
                  
        end
    end
end

P                 =  graph(connectivity,'upper');
bins              =  conncomp(P);
binnodes = accumarray(bins', 1:numel(bins), [], @(v) {sort(v')});
% if numel(binnodes) == 1 
%     fprintf('All landmarks are connected.');
% else
%     fprintf('number of Regions = %d\n\n', numel(binnodes));
% end
% fprintf('\n');

bin_N = numel(binnodes);

connectivity      =  connectivity + connectivity' ;
K                 =  length(chart);

connectivity(1:K+1:end)=ones(K,1);   % denote each node is connected with itself. The distance is set as 1 just for convention.



end


