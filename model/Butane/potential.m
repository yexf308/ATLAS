function [V] = potential(x1,y1,y3,x4,y4,z4,parameter)
        l = parameter.l;
        k2 = parameter.k2;
        k3 = parameter.k3; 
        theta = parameter.theta;
        c1   = parameter.c1;
        c2   = parameter.c2;
        c3   = parameter.c3;
        
        b1 = sqrt(x1.^2 + y1.^2) - l;
        b2 = abs(y3) - l; 
        b3 = sqrt(x4.^2 + z4.^2 + (y3 - y4).^2) - l;
        
        theta1 = theta - acos( y1.*sign(y3)./sqrt(x1.^2 + y1.^2) );
        theta2 = theta - acos( (y3 - y4).*sign(y3)./sqrt( x4.^2 +  z4.^2 + (y3 - y4).^2 ));
        
        c1_term = c1 * x4.*sign(x1)./(x4.^2 + z4.^2).^(1/2);
        c2_term = c2 * x4.^2./(x4.^2 + z4.^2);
        c3_term = c3 * x4.^3 .* sign(x1)./(x4.^2 + z4.^2).^(3/2);
        
        V = c1_term + c2_term + c3_term + 1/2 * k2 * (b1.^2 + b2.^2 + b3.^2) + 1/2 * k3 * (theta1.^2 + theta2.^2);
        

end

