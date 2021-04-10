function [Rx] = Rtrue(X,R)
         x = X(1);
         y = X(2);
         z = X(3);
         theta           = atan2( sqrt(x.^2+y.^2),z);
         phi             = atan2(  y, x);
         Rx = R(theta, phi);
         
end