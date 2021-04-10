function [diffusion] = get_diffusion(X,parameter)

sigma   = parameter.sigma;

diffusion = sigma;
%diffusion = sigma + sqrt(lambda) *(y - tau*sin(omega*x) )^2 ;

end