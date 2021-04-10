function [b] = btrue(X,b_ito, b_stra)
         x = X(1);
         y = X(2);
         z = X(3);
         theta           = atan2( sqrt(x.^2+y.^2),z);
         phi             = atan2(  y, x);
         b               = b_ito(theta,phi) + b_stra(theta,phi);  
         
         
end

