function c = correlation_matern ( Nx, x, rho, nu )
% Computes Matern correlation with smoothness nu and scale rho for
% separations x, following definition on wikipedia:
%   r = ( sqrt(2*nu) * x / rho ) 
%   c =  r.^nu .* besselk( nu, r ) * ( 2^(1-nu) / gamma(nu) )
%
% Note this has special cases:
% nu =  1/2 ,  c = exp( - r / rho )   
% nu ? infty,  c = exp( - 0.5 * r.^2 / rho^2 )   
%
% To reproduce previous implementation, new rho should be reduced by factor
% of 1/sqrt(2) relative to previous rho0.

if min( x ) < 0
    disp([' in correlation_matern, min(x) is < 0,  = ' num2str(min(x)) ])
return
end

scaledr = ( sqrt(2*nu) * x / rho ) ;
c =  ( scaledr.^nu .* besselk( nu, scaledr ) * ( 2^(1-nu) / gamma(nu) ) )' ;
ii = find( x == 0 ); c(ii) = 1. ; % fix separation = 0: besselk(nu,0) returns NaN

return
