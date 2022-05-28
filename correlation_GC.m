function c = correlation_GC(Nx,x,cutoff,dummy)
% Gaspari-Cohn correlation fn. (periodic properties??)
% From(4.10) of Gaspari and Cohn (QJRMS, 125, 723-757).
% Implementation here copied from compact_correlation.m
% Inputs: Nx = length(x)
%         x  = separations for evaluation of correlation
%              ... in principle, these could be arbitrary, but, for 
%              convenience in periodic domain, assume evenly spaced with
%              dx = domain_length / Nx.
%     cutoff = length scale for periodic, squared exponential correlation
%      dummy = dummy argument, for calls from makeCorrelationMatrix.m
%
% A QUESTION: have notes from long ago (test_2scales_....m) stating that GC
% is not valid correlation fn on periodic domain. True?  (If so, what about
% all the applications to L96?? Source of instability in serial update?)

if min( x ) < 0
    disp([' in correlation_GC, min(x) is < 0,  = ' num2str(min(x)) ])
return
end

c = zeros( size( x ) )';  % ' , because want output to be col vector

indices = find( x >= cutoff / 2 & x < cutoff );
r = 2 * x(indices) / cutoff ;
c(indices) =  r.^5 / 12 - r.^4 / 2 + r.^3 * 5/8 + r.^2 * 5/3   ...
-r .* 5 + 4 - 2 ./ (3 * r);

indices = find( x < cutoff / 2 );
r = 2 * x(indices) / cutoff ;
c(indices) = -r.^5 / 4 + r.^4 / 2 + r.^3 * 5/8 - r.^2 * 5/3 + 1;

return
