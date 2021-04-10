Lambda_true        = Exact_parameter.Lambda_true;
b_true             = Exact_parameter.b_true;
R_true             = Exact_parameter.R_true;
D                  = Exact_parameter.D;

K                  = length(chart);

absolute_error     = zeros(4, K);
relative_error     = zeros(2, K);

for i = 1:K
    X0                    = chart{i}.X_int;
    U                     = chart{i}.U;
    Lambda_hat            = chart{i}.Lambda;
    Lambda_X0_true        = Lambda_true(X0);
    [U_true,Sigma_square] = svd( Lambda_X0_true );
   
    
    U_true                = U_true(:,1:d);
        
    b                     = chart{i}.b;    
    b_X0_true             = b_true(X0);
    
    absolute_error(1, i)  = norm( Lambda_hat - Lambda_X0_true );  % Fro-norm of H is sqrt(sum(sigma.^2))
    absolute_error(2, i)  = norm( b_X0_true - b );
    absolute_error(3, i)  = norm( norm(X0) - R_true(X0) );
    absolute_error(4, i)  = subspace(U,U_true);   % angle between two subspace
    
    relative_error(1, i)  = absolute_error(1, i) /norm(Lambda_X0_true);
    relative_error(2, i)  = absolute_error(2, i)/norm( b_X0_true );
    
end


%[m1, i1]                                  = max(absolute_error(1,:));
%[m2, i2]                                  = max(absolute_error(2,:));
%[m3, i3]                                  = max(absolute_error(3,:));
%[m4, i4]                                  = max(absolute_error(4,:));
%disp(['maximum absolute error of Fro-norm of H is ', num2str(m1), ' at chart No.', num2str(i1)])
%disp(['maximum absolute error of b is ', num2str(m2), ' at chart No.', num2str(i2)])
%disp(['maximum absolute error of landmark is ', num2str(m3), ' at chart No.', num2str(i3)])
%disp(['maximum error of angle of subspace of U is ', num2str(m4), ' at chart No.', num2str(i4)])
%[b1, i4]                                  = max(norm_b);
%[H1, i5]                                  = max(norm_H);
%disp(['maximum norm of b is ', num2str(b1), ' at chart No.', num2str(i4)])
%disp(['maximum norm of H is ', num2str(H1), ' at chart No.', num2str(i5)])


