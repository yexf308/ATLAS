function [nonlin_trans] = get_nonlin_trans(X,parameter)
%         D     = length(X);
%         theta = X(1);
%         r = X(2);
%         nonlin_trans = zeros(D,1);
%         t = parameter.t;
%         Tran = parameter.Tran;
%         
%                
%         
%         nonlin_trans(1) = r.*cos(theta+r+t);
%         nonlin_trans(2) = r.*sin(theta+r+t);
%         nonlin_trans(3:end) = Tran*X(2:end); 
t     = parameter.t;
Tran  = parameter.Tran;

if length(size(X)) == 2
     % X is D*tn matrix

    [D,n] =size(X); 
    theta = X(1,:);
    r     = X(2,:);
    nonlin_trans = zeros(D,n);

    nonlin_trans(1,:) = r .* cos(theta + r + t);
    nonlin_trans(2,:) = r .* sin(theta + r + t);
    nonlin_trans(3:end,:) = Tran * X(2:end, :);
else
    %X is N*D*tn tensor
    
%     [N, D, tn] = size(X);
%     theta      = X(:,1,:);
%     r          = X(:,2,:);
%     nonlin_trans  = zeros(N,D,tn);
%     
%     nonlin_trans(:,1,:) = r .* cos(theta + r + t);
%     nonlin_trans(:,2,:) = r .* sin(theta + r + t);
%     nonlin_trans(:,3:end,:) = permute( reshape( Tran * reshape( permute(X(:,2:end,:), [2,1,3]), D-1, []) , D-2,N, tn), [2,1,3]);
end

end

