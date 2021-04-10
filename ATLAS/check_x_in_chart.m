function [result] = check_x_in_chart(x, chart, t0, r, chi_p,threshold)
% This function checks whether x is in the chart of time r*t0. Usually,
% 0<r<4, 0<d(x,z_L)<2sqrt(t0)

% chi_p is 95% statistics, in 2D, it is 5.991
if ~isrow(x)
    x=x';  % make sure x is a row
end

sigma       = chart.sigma;
X_int       = chart.X_int;
WU          = chart.WU;

Center      = X_int; 
x_coeff     = (x - Center)*WU; 
result      = 0;
if  sum(x_coeff.^2./sigma.^2) <= r^2 * t0 * chi_p ...
 && (norm(x-X_int)<threshold )  %this is  an very ad-hoc threshold. which is determined by the thickness of fast mode.
                                 % this is necessary when building the
                                 % connectivity of charts. some charts that are far away in
                                 % Euclidean distance may be very close in
                                 % diffusion distance. 
    result=1;  
    
end

end

