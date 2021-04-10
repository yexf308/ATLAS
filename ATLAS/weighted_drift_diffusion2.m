function [ X_proj, b, Lambda_hat, H_hat, T, nearest, neigh ] = weighted_drift_diffusion2( X0, chart, neigh, nearest, connectivity_indices, t0, chi_p, D,d, threshold, option,mode )

index           = 1;
while index
        L_n             = length(neigh);
        T               = zeros(1, L_n);
        for k = 1:L_n
            cur_chart        = chart{neigh(k)};
            Center           = cur_chart.X_int;
            sigma            = cur_chart.sigma;
            if option == 1    
                 WU               = cur_chart.WU; % this is fast mode projection
                 x_coeff          = ( X0 - Center )*WU;

             elseif option == 2
                 U                = cur_chart.U;% this is orthorgnal projection
                 x_coeff          = ( X0 - Center) * U;
             end
            t                = sum( x_coeff.^2 ./ sigma.^2 )/(chi_p);
            T(k)             = min( sqrt(t/t0), 10 );
        end
        [~, i] = min(T);
        if neigh(i) == nearest
            index = 0;
        else
            nearest = neigh(i);
            neigh   = connectivity_indices{nearest}; 
        end
end
    

L_n             = length(neigh);
T               = zeros(1, L_n);
Lambda          = zeros(D, D);
b               = zeros(D, 1);
X               = zeros(L_n,D);
weight          = 0;
X_proj          = X0;

%% Project X_curr on the estimated manifold and estimate the effetive drift and diffusion
for k = 1:L_n
    cur_chart        = chart{neigh(k)};
    Center           = cur_chart.X_int;
    if mode == 1 && norm(X0-Center) > threshold(1)                                       %  if mode = explore and out of chart
         T(k)        = 10;
         X(k, :)     = zeros(1, D);
            % X(k,:) remains zeros, weight contribution for this is also zero.    
    else
    
         U                = cur_chart.U;
         Ut               = U';
         sigma            = cur_chart.sigma;
         

         if option == 1 
             WU               = cur_chart.WU;                        % this is fast mode projection
             x_coeff          = ( X0 - Center )*WU;
            
         elseif option == 2                                           % this is orthorgnal projection
             x_coeff     = ( X0 - Center) *U;
         end
        
         t           = sum( x_coeff.^2 ./ sigma.^2 )/(chi_p);
         T(k)        = min( sqrt(t/t0), 10 );
         expmTk      = exp( -T(k) );
         X(k,:)      = ( x_coeff * Ut + Center )  * expmTk;

         weight      = weight + expmTk;
         b           = b      + cur_chart.b * expmTk;
         Lambda      = Lambda + cur_chart.Lambda * expmTk;

    end
end


if weight~=0
    X_proj              = sum(X,1)/weight;
    b                   = b/weight;
    Lambda              = Lambda/weight;
end

[UL, SL, ~]         = svd(Lambda);
UL                  = UL(:,1:d);
SL                  = SL(1:d,1:d);
Lambda_hat          = UL * SL * UL';
H_hat               = UL * sqrt( SL );

% Likely useful when D will be much larger; the LambdaUd can be computed once for all in each chart, and interpolated
% [U2,S2,~]               = svds(Lambda,d,'largest','SubspaceDimension',d+2,'Tolerance',1e-3,'LeftStartVector',?LambdaUd?,'MaxIterations',10);%,'Display',true);
% Lambda_hat2             = U2(:,1:d) * S2(1:d, 1:d) * U2(:, 1:d)';
% H_hat2                  = U2(:,1:d) * sqrt( S2(1:d, 1:d) );
%norm(Lambda_hat-Lambda_hat2)/norm(Lambda_hat)
% min(norm(H_hat-H_hat2),norm(H_hat+H_hat2))

end

