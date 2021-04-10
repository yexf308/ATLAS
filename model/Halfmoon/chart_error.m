K = length(chart);

angle          = zeros(1,K);
manifold_error = zeros(1,K);
a1             = parameter.a1;
a2             = parameter.a2;
a3             = parameter.a3;
a4             = parameter.a4;
relative_error = zeros(4,K);
theta_s        = atan(b2^2/(2*b1));
r_true         = exp(-b2^2/(4*b1))*sqrt(1+(b2^2/(2*b1))^2);

for i = 1:K
    X_int             = chart{i}.X_int;
    z1                = X_int(1);
    z2                = X_int(2);
    U                 = chart{i}.U;
    H_hat             = chart{i}.U*chart{i}.sigma;
    
    theta_curr        = atan2(z2, z1);
    theta             = theta_curr - theta_s;
    
    
    manifold_error(i) = sqrt(   ( r_true - sqrt(z1^2 + z2^2))^2 + sum((X_int(3:end)-1).^2) ); 
    b_stra            = a1 + a2 * sin(2*theta) + a3* cos(theta);
    b1                = -z2 .*b_stra - z1.*a4^2/2;
    b2                = z1 .*b_stra - z2.*a4^2/2;
    b_true            = [b1; b2; zeros(D-2,1)];
    H_true            = [ -a4*z2; a4*z1;zeros(D-2,1)]; 
    Lambda_true       = H_true * H_true';
    Lambda_hat        = H_hat * H_hat';
    
    angle(i)            = subspace(H_true, H_hat);
    sigma               = chart{i}.sigma;
    b                   = chart{i}.b;
    relative_error(1,i) = norm(b(1:2) - b_true(1:2))/norm(b_true(1:2)); % first two 
    relative_error(2,i) = norm(b - b_true)/norm(b_true); %total
    relative_error(3,i) = norm( Lambda_true(1:2,1:2) - Lambda_hat(1:2,1:2))/norm(Lambda_true(1:2,1:2)); %first two 
    relative_error(4,i) = norm( Lambda_true - Lambda_hat)/norm(Lambda_true); %total
    
end
