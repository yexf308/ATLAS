function [Lambda] = Lambdatrue(X,H_true)
         x = X(1);
         y = X(2);
         z = X(3);
         theta           = atan2( sqrt(x.^2+y.^2),z);
         phi            = atan2(  y, x);
         
         Lambda=H_true(theta, phi) * H_true(theta, phi)';
         
end