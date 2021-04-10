function [diffusion] = get_diffusion(X,parameter)
  
        D         = length(X);
        diffusion = zeros(D,D);
        a4 = parameter.a4;
        b2 = parameter.b2; 
        b4 = parameter.b4;
        
       diffusion(1,1) = a4;
       
       diffusion(D+2:D+1:end) = ones(1,D-1) * b4;
       diffusion(2,2) = b2;
    
   
end
