function [chart] = Learning_Slow_Manifold(data, Cov_store, Mean_store, D, d , modify )
%% Set the parameter 

dt              = data.dt;
UpperBound      = data.UpperBound;
LowerBound      = data.LowerBound;

s_u             = round(UpperBound/dt);
s_l             = round(LowerBound/dt)+1;

s_span          = s_l:s_u;
t_span          = s_l*dt:dt:s_u*dt;
t_span          = t_span - mean(t_span);

Cov_span        = Cov_store(s_span);
Mean_span       = Mean_store(s_span);

sum_Cov         = zeros(D,D);
sum_Mean        = zeros(1,D);

M               = length(s_span);

%% Calculate the diffusion matrix Lambda and drift b explicitly
for i = 1:M
    sum_Cov      = sum_Cov  + Cov_span{i};
    sum_Mean     = sum_Mean + Mean_span{i};
end

temp_Cov          = zeros(D,D);
temp_Mean         = zeros(1,D);

for i =1:M
    Cov_span{i}        = Cov_span{i} - sum_Cov/M; 
    temp_Cov           = temp_Cov  + Cov_span{i}*t_span(i);
    
    Mean_span{i}       = Mean_span{i} - sum_Mean/M;
    temp_Mean          = temp_Mean + Mean_span{i}*t_span(i); 
end

%b_hat = \sum(m_bar(t)*t_bar)/|t_bar|^2;
%Lambda_hat =  \sum(Lambda_bar(t)*t_bar)/|t_bar|^2;

Lambda_full                 = temp_Cov ./norm(t_span)^2;
b_full                      = temp_Mean./norm(t_span)^2;

if modify == 0
    x0                          = sum_Mean/M - b_full * (UpperBound + LowerBound)/2;   % Use intercept as Landmark of the chart
else
    x0                          = sum_Mean/M - b_full * (UpperBound - LowerBound)/2;   % Use tau_min time point as Landmark of the chart
end
Gamma_full                  = sum_Cov/M  - Lambda_full * (UpperBound + LowerBound)/2;
%% Projection to d dimension

[U, S, ~]              = svd(Lambda_full);
Lambda_hat             = U(:,1:d)*S(1:d, 1:d)*U(:,1:d)';
sigma                  = sqrt(diag(S(1:d, 1:d)));
sigma                  = sigma'; %make sure sigma is a row vector
U                      = U(:,1:d);


[V,VS, ~]              = svd(Gamma_full);
diag_VS                = diag(VS);
index                  = find(diag_VS>(diag_VS(1)/3));
V_index                = [];
for i = 1: length(index)
    if subspace(U, V(:,index(i)))/pi >0.2  %this need to be changed
        V_index        = [V_index, V(:,index(i))];
    end
end


E                      = [U, V_index];
[Q,R]=qr(E);
WU    = Q*pinv(R*R')*Q'*U;




%% Store data into the chart
 chart = struct(                  ...
      'LearningTime', NaN,        ...
      'U',            U,          ...
      'V',            V_index,    ...
      'sigma',        sigma,      ...
      'sigma_fast',   diag(VS),   ...
      'Lambda',       Lambda_hat, ...
      'b',            b_full',    ...   % b is a column vector
      'X_int',        x0,         ...   % X_int is a row vector
      'WU',           WU          ...
          );

end

