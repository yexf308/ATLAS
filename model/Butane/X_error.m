function [manifold_error,b_error, b_abs_error, sigma_error,sigma_abs_error,angle ] = X_error(chart, parameter)
    X     = chart.X_int;
    sigma = chart.sigma;
    U     = chart.U;
    b     = chart.b;
    H     = sigma*U;
    x1    = X(1);
    y1    = X(2);
    y3    = X(3);
    x4    = X(4);
    y4    = X(5);
    z4    = X(6);   
    l     = parameter.l;
    theta = parameter.theta;
    c1    = parameter.c1;
    c2    = parameter.c2;
    c3    = parameter.c3;
    manifold_error = sqrt((X(1) + l*sin(theta) ) .^2 + ( X(2) - l*cos(theta) ).^2 + (X(3) - l ).^2 + (X(5) +l*cos(theta) -l ).^2 ...
          + (sqrt(X(4).^2 + X(6).^2) - l*sin(theta)).^2 ); 
      % note this is not the 2-norm error ) 
      % true invariant manifold is [l*sin(theta), - l*cos(theta), -l, x,
      % l*cos(theta)-l, z], with x^2+z^2 = l^2sin(theta)^2.
    sigma_true =   parameter.sigma/(l*sin(theta)); 
    
    H_true  = [0  0  0 -z4.*sigma_true 0 x4.*sigma_true]';
    sigma_abs_error = norm(norm(H_true,'fro') - abs(sigma));
    sigma_error  =  sigma_abs_error/norm(H_true,'fro');
    angle = subspace(H_true, H);
    
    c_term =  z4.*   (  (c1 .*(x4.^2 + z4.^2 )  + 3 * c3 *x4.^2 ) + 2* c2*x4.*sqrt(x4.^2 + z4.^2) )./(x4.^2 + z4.^2).^(5/2);
   %c_term =  z4.*   (  (c1 .*(l*sin(theta))^2  + 3 * c3 *x4.^2 ) + 2* c2*x4.*(l*sin(theta)) )./(l*sin(theta)).^5;
   
  b4    =  - c_term * z4 - sigma_true^2/2*x4;
  b6    =    c_term * x4 - sigma_true^2/2*z4;
  
  
  %  c_term = -c1*z4./sqrt(x4.^2 + z4.^2) - 2*c2 * x4.*z4./(x4.^2 + z4.^2) - 3*c3 * x4.^2.*z4./(x4.^2 + z4.^2).^(3/2);
  %  b4    =   c_term * z4 - sigma_true^2/2*x4;
  % b6    =   - c_term * x4 - sigma_true^2/2*z4;
   
  b_true       = [ 0 0 0 b4 0 b6]';
  b_abs_error  =  norm(b-b_true);
  b_error      = norm(b-b_true)/norm(b_true);
 % b_abs_error  =  norm(b([4,6])-b_true([4,6]));
 % b_error      = norm(b([4,6])-b_true([4,6]))/norm(b_true);
  
end

