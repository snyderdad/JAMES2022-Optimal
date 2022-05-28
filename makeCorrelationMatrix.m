function [ C ] = makeCorrelationMatrix_new( x, x0, nu, correlation )
% Construct a correlation matrix
%
% Inputs:
%   x = locations as a row vector
%   x0 = length scale
%   nu = Matern parameter (in Matern case ... can use other fns, params)
%   @CORRELATION, a handle for a correlation function.
%
% Output:
%   C = correlation matrix for input locations x
%

Nx = length(x) ;
% matrix of indices s.t. x(Ilocs(i,j)) is separation between ith, jth locations
Ilocs = abs( repmat( 1:Nx, [Nx,1]) - repmat( (1:Nx)', [1,Nx]) ) + 1 ;

C = correlation( Nx, x(Ilocs), x0, nu ) ;
                % at this point, code for correlation fns has evolved to where
                % 1st argument irrelevant; but a pain to change this interface 

end
