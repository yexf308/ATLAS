function [drift] = get_drift(X,parameter)
        D     = length(X);
        theta = X(1);
        r  = X(2);
        r2 = X(3:end);
        drift = zeros(D,1);
        a1 = parameter.a1;
        a2 = parameter.a2;
        a3 = parameter.a3;
        b1 = parameter.b1;
        b3 = parameter.b3;
        
        
        
        drift(1)      = a1+a2*sin(2*theta)+a3*cos(theta);
        drift(2)      = b1.*(1-r);
        drift(3:end)  = -b3.*r2;
        
        
        

end
