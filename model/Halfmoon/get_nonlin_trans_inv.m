function [nonlin_trans_inv] = get_nonlin_trans_inv(X,parameter)
%         x = X(1);
%         y = X(2);
%         D = length(X);
%         nonlin_trans_inv = zeros(D,1);
%         t = parameter.t;
%         Tran_inv = parameter.Tran_inv;
%         angle =  mod(atan2(y,x),2*pi); 
%         r = sqrt(x.^2+y.^2);
%         if angle - r - t>pi
%              nonlin_trans_inv(1) = angle - r - t-2*pi;
%         else
%              nonlin_trans_inv(1) = angle - r - t;
%         end
%         
%         nonlin_trans_inv(2) = r;
%         nonlin_trans_inv(3:end) = Tran_inv * [r;X(3:end)];
    t               = parameter.t;
    Tran_inv        = parameter.Tran_inv;
    
if length(size(X)) == 2
    x               = X(1,:);
    y               = X(2,:);
    [D, n]          = size(X);
    nonlin_trans_inv = zeros(D,n);
    angle           =  mod(atan2(y,x),2*pi); 
    r               = sqrt(x.^2+y.^2);

    nonlin_trans_inv(1,:) = angle - r - t;
    index = find((angle-r-t)>pi);
    nonlin_trans_inv(1,index) = angle(index) - r(index) - t-2*pi;
    nonlin_trans_inv(2,:) = r;
    nonlin_trans_inv(3:end,:) = Tran_inv * [r;X(3:end,:)];
end
     %X is N*D*tn tensor
     
%     [N, D, tn] = size(X);
%     x               = X(:,1,:);
%     y               = X(:,2,:);
%     nonlin_trans_inv = zeros(N,D,tn);
%     angle            = mod(atan2(y,x),2*pi); 
%     r                = sqrt(x.^2+y.^2);
%     
%     M                = angle - r - t;
%     M_col            = reshape(M, N*tn,1);
%     index            = find(M_col>pi);
%     M_col(index)     = M_col(index) - 2*pi;
%     M                = reshape(M_col, N,1,tn);
%     
%     nonlin_trans_inv(:,1,:) = M;
%     nonlin_trans_inv(:,2,:) = r;
%     nonlin_trans_inv(:,3:end,:) = permute( reshape( Tran_inv * reshape( permute([r, X(:,2:end,:)], [2,1,3]), D-1, []) , D-2,N, tn), [2,1,3]);
%     
    
    
end