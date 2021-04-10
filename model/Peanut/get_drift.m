function [drift] = get_drift( X,epsilon,a1,a2,c1,c2,c3,c4,c5,c6 )

x  = X(1);
y  = X(2);
z  = X(3);

r           = sqrt(x.^2+y.^2+z.^2);
r_sq        = r.^2;
rsin        = sqrt(x.^2 + y.^2);
rsin_sq     = rsin.^2;  %x.^2+y.^2


J           = [ x./r, x.*z./rsin, -y;
                y./r, y.*z./rsin,  x;
                z./r, -rsin,       0 ];
            
b           = [ - c1/epsilon * ( 1 - sqrt( a1 + a2 * z.^2./r_sq ) ./ r ) ;
                  c3 * ( 4*z.^3./r.^3-3*z./r )./rsin                     ;   
                  c5 * ( y.*z./(rsin.* r_sq) + x./r_sq )                ];
                 
%Ito         = 1/2* ([ -x; -y; -z]* c4^2 * rsin_sq./(r_sq.^2) + [-x; -y; 0] * c6^2./r_sq);
d1          = c4^2 * rsin_sq./(r_sq.^2);
d2          = c6^2./r_sq;
sumd1d2     = d1+d2;
Ito         = 1/2 * [-x * sumd1d2 ; -y* sumd1d2; -z *d1  ];

drift       = J*b + Ito;



% x  = X(1);
% y  = X(2);
% z  = X(3);
% 
% r           = sqrt(x.^2+y.^2+z.^2);
% r_sq        = r.^2;
% rsin        = sqrt(x.^2 + y.^2);
% rsin_sq     = rsin.^2;
% 
% 
% J           = [ x./r, x.*z./rsin, -y;
%                 y./r, y.*z./rsin,  x;
%                 z./r, -rsin,       0 ];
% 
% b           = [ - c1/epsilon * ( 1 - sqrt( a1 + a2 * z.^2./r_sq ) ./ r ) ;
%                   c3 * ( 4*z.^3./r.^3-3*z./r )./rsin                     ;
%                   c5 * ( y.*z./(rsin.* r_sq) + x./r_sq )                ];
% 
% Ito         = 1/2* ([ -x; -y; -z]* c4^2 * rsin_sq./(r_sq.^2) + [-x; -y; 0] * c6^2./r_sq);
% 
% drift       = J*b + Ito;

end