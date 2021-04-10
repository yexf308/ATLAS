function [drift] = get_drift(X,parameter)
        x1 = X(1);
        y1 = X(2);
        y3 = X(3);
        x4 = X(4);
        y4 = X(5);
        z4 = X(6);
        drift = zeros(6,1);
        l = parameter.l;
        k2 = parameter.k2;
        k3 = parameter.k3; 
        theta = parameter.theta;
        c1   = parameter.c1;
        c2   = parameter.c2;
        c3   = parameter.c3;
        
        b1 = 1 - l./sqrt(x1.^2 + y1.^2);
        b3 = 1 - l./sqrt(x4.^2 + z4.^2 + (y3 - y4).^2);
        
        theta1 = theta - acos( y1.*sign(y3)./sqrt(x1.^2 + y1.^2) );
        theta2 = theta - acos( (y3 - y4).*sign(y3)./sqrt( x4.^2 +  z4.^2 + (y3 - y4).^2 ));
        
        %c_term = - z4.*    ( sign(x1).*(c1 .*(x4.^2 + z4.^2 )  + 3 * c3 *x4.^2 ) + 2* c2*x4.*sqrt(x4.^2 + z4.^2) )./(x4.^2 + z4.^2).^(5/2);
        
        c_term =  z4.*    (  (c1 .*(x4.^2 + z4.^2 )  + 3 * c3 *x4.^2 ) + 2* c2*x4.*sqrt(x4.^2 + z4.^2) )./(x4.^2 + z4.^2).^(5/2);
        
        drift(1) = k2 * x1 .*b1        -  k3 * y1.*sign(x1.*y3)./(x1.^2 + y1.^2 ).* theta1;
        
        drift(2) = k2 * y1 .*b1        +  k3 * abs(x1).*sign(y3)./(x1.^2 + y1.^2 ).* theta1;
        
        drift(3) = k2 * ((y3 - y4).* b3 + (-l + abs(y3)).*sign(y3) )+  k3 * sqrt(x4.^2 +z4.^2).* sign(y3)./(x4.^2 + z4.^2 + (y3 - y4).^2).*theta2;
        
        drift(4) = k2 * x4 .* b3       -  k3 * x4.*(y3-y4)./sqrt(x4.^2 +z4.^2).*sign(y3)./(x4.^2 + z4.^2 + (y3 - y4).^2).*theta2  ...
                    +z4.*c_term;
                
        drift(5) = -k2 * (y3 - y4).* b3  -  k3 * sqrt(x4.^2 +z4.^2).* sign(y3)./(x4.^2 + z4.^2 + (y3 - y4).^2).*theta2;
        
        drift(6) = k2 * z4 .* b3       -  k3 * z4.*(y3-y4)./sqrt(x4.^2 +z4.^2).*sign(y3)./(x4.^2 + z4.^2 + (y3 - y4).^2).*theta2  ...
                   -x4.*c_term;
        
        drift = -drift;

end

