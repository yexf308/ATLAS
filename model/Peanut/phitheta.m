function [index] = phitheta(X_curr, phi_zero, theta_zero)
       theta_curr           = atan2( sqrt( X_curr(1)^2 + X_curr(2)^2 ), X_curr(3) );
       phi_curr             = mod(atan2(X_curr(2), X_curr(1)),2*pi);
       
       
       [~,i]=min(abs(phi_curr-phi_zero));
       if theta_curr>theta_zero(i)
           index = -1;
       else
           index = 1;
       end

end

