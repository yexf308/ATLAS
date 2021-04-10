% K                   = length(chart);
% manifold_error      = zeros(1,K);
% b_error             = zeros(1,K);
% b_abs_error         = zeros(1,K);
% norm_H_Fro_error    = zeros(1,K);
% norm_H_Fro_abs_error= zeros(1,K);
% sigma_fast          = zeros(6,K);
% fast_slow_angle      = zeros(1,K);
% U_angle              = zeros(1,K);
% norm_b             = zeros(1, K);
% norm_H             = zeros(1, K);
% R_error            = zeros(1, K);
% phi_store          = zeros(1, K);
% X_int_store        = zeros(K, 6);
% for i = 1:K
%     [manifold_error(i),b_error(i), b_abs_error(i), norm_H_Fro_error(i),norm_H_Fro_abs_error(i), U_angle(i)] = X_error(chart{i}, parameter);
%     sigma_fast(:,i)       =chart{i}.sigma_fast;
%     X_int                 = chart{i}.X_int;
%     X_int_store(i,:)      = X_int;
%     phi_store(i)          = mod(atan2(X_int(6), X_int(4)),2*pi);
%     R_error(i)            = sqrt(X_int(4)^2 + X_int(6)^2) - l*sin(theta);
%     b                     = chart{i}.b;
%     norm_b(i)             = norm(b);
%     U                     = chart{i}.U;
%     V                     = chart{i}.V;
%     sigma                 = chart{i}.sigma;
%     H                     = U*diag(sigma);
%     norm_H(i)             = norm(H);
%     fast_slow_angle(i)    = subspace(V, U)/pi;
% end
% 
% absolute_error      = zeros(4, K);
% absolute_error(1,:) = norm_H_Fro_abs_error;  % Fro-norm of H is abs(sigma)
% absolute_error(2,:) = b_abs_error;            
% absolute_error(3,:) = manifold_error;
% absolute_error(4,:) = U_angle;
% 
% relative_error      = zeros(2,K);
% relative_error(1,:) = norm_H_Fro_error;
% relative_error(2,:) = b_error;
% 

%[m1, i1]                                  = max(absolute_error(1,:));
%[m2, i2]                                  = max(absolute_error(2,:));
%[m3, i3]                                  = max(absolute_error(3,:));

%[m4, i4]                                  = max(absolute_error(4,:));
%disp(['maximum absolute error of Lambda is ', num2str(m1), ' at chart No.', num2str(i1)])
%disp(['maximum absolute error of b is ', num2str(m2), ' at chart No.', num2str(i2)])
%disp(['maximum absolute error of landmark is ', num2str(m3), ' at chart No.', num2str(i3)])
%disp(['maximum absolute error of tangent line is ', num2str(m4), ' at chart No.', num2str(i4)])
%[b1, i4]                                  = max(norm_b);
%[H1, i5]                                  = max(norm_H);
%disp(['maximum norm of b is ', num2str(b1), ' at chart No.', num2str(i4)])
%disp(['maximum norm of H is ', num2str(H1), ' at chart No.', num2str(i5)])

% angle between top 2 eigenvector contribution from fast mode and the slow
% mode is almost pi/2. Also the top 2 eigenvalues are the same (0.0066). 
% But the the third eigenvector of the diffusion
% matrix is not orthgonal (many times almost parallel) with the slow mode. 

K                   = length(chart);
manifold_error      = zeros(1,K);
angle               = zeros(1,K);
relative_error      = zeros(4,K);
l                   = parameter.l;
theta               = parameter.theta;
c1                  = parameter.c1;
c2                  = parameter.c2;
c3                  = parameter.c3;


for i = 1:K
    X_int  = chart{i}.X_int;
    
    manifold_error(i) = sqrt((X_int(1) + l*sin(theta) ) .^2 + ( X_int(2) - l*cos(theta) ).^2 + (X_int(3) - l ).^2 + (X_int(5) +l*cos(theta) -l ).^2 ...
          + (sqrt(X_int(4).^2 + X_int(6).^2) - l*sin(theta)).^2 ); 
    
    U                 = chart{i}.U;
    H_hat             = chart{i}.U*chart{i}.sigma;
    
    x4         = X_int(4);
    z4         = X_int(6);
    sigma_true =   parameter.sigma/(l*sin(theta));    
    H_true     = [0  0  0 -z4.*sigma_true 0 x4.*sigma_true]';
    
    angle(i)   = subspace(H_true, H_hat);
    
    c_term     =  z4.*   (  (c1 .*(x4.^2 + z4.^2 )  + 3 * c3 *x4.^2 ) + 2* c2*x4.*sqrt(x4.^2 + z4.^2) )./(x4.^2 + z4.^2).^(5/2);
    b4         =  - c_term * z4 - sigma_true^2/2*x4;
    b6         =    c_term * x4 - sigma_true^2/2*z4;
    
    b_true       = [ 0 0 0 b4 0 b6]';
    b            = chart{i}.b;
    
    Lambda_true = H_true * H_true';
    Lambda_hat  = H_hat  * H_hat';
    
    relative_error(1,i) =  norm(b - b_true)/norm(b_true);
    relative_error(2,i) =  norm(Lambda_true - Lambda_hat)/norm(Lambda_true);
    relative_error(3,i) =  norm(b([4,6]) - b_true([4,6]))/norm(b_true([4,6]));
    relative_error(4,i) =  norm(Lambda_true([4,6],[4,6]) - Lambda_hat([4,6],[4,6]))/norm(Lambda_true([4,6],[4,6]));
end