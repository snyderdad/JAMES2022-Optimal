% Compute MSE, averaged over Nreps realizations, for std 2-scale example
% after "obvious" transformation (which is the optimal, since R = H = I,
% here) and complete localization (i.e. assume all cross-covs are zero).
% For manuscript, where it will be compared to MSE from usual localization.
% 
% Code cannibalized from test_updateWithOptlTransf.m, 28 June 2019

% std 2-scale example: 
%   P = sigma_l^2 * 1 + sigma_s^2 * I ; H = R = I

Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^4 ; obs_are = 'EquallySpaced' ;
sigma2_s = 0.25 ; sigma2_l = 1 ; r2 = 2 ;
fig_number = 106 ;

Nx = 25 ; Ny = Nx ; Ne = 10 ; Nreps = 10^3 ; obs_are = 'EquallySpaced' ;
sigma2_s = 1 ; sigma2_l = 1 ; r2 = 2 ;
fig_number = 107;

Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^3 ; obs_are = 'EquallySpaced' ;
sigma2_s = 1 ; sigma2_l = 1 ; r2 = 2 ;
fig_number = 108;

Nx = 100 ; Ny = Nx ; Ne = 4 ; Nreps = 10^3 ; obs_are = 'EquallySpaced' ;
sigma2_s = 1 ; sigma2_l = 1 ; r2 = 2 ;
fig_number = 109;

Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^3 ; obs_are = 'EquallySpaced' ;
sigma2_s = 1 ; sigma2_l = 1 ; r2 = 4 ;
fig_number = 110;

Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^3 ; obs_are = 'EquallySpaced' ;
sigma2_s = 2 ; sigma2_l = 1 ; r2 = 2 ;
fig_number = 111;

Nx = 500 ; Ny = Nx ; Ne = 4 ; Nreps = 10^4 ; obs_are = 'EquallySpaced' ;
sigma2_s = 1 ; sigma2_l = 1 ; r2 = 2 ;
fig_number = 113;

Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^4 ; obs_are = 'EquallySpaced' ;
sigma2_s = 1 ; sigma2_l = 1 ; r2 = 2 ;
fig_number = 13;

Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^4 ; obs_are = 'EquallySpaced' ;
sigma2_s = 1 ; sigma2_l = 0.5 ; r2 = 2*(sigma2_s + sigma2_l) ;
fig_number = 14;

%Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^4 ; obs_are = 'EquallySpaced' ;
%sigma2_s = 1 ; sigma2_l = 0.5 ; r2 = 0.5*(sigma2_s + sigma2_l) ;
%fig_number = 15;

Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^5 ; obs_are = 'EquallySpaced' ;
sigma2_s = 0.5 ; sigma2_l = 1.0 ; r2 = 0.5*(sigma2_s + sigma2_l) ;
fig_number = 16;

Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^5 ; obs_are = 'EquallySpaced' ;
sigma2_s = 0.5 ; sigma2_l = 1.0 ; r2 = 0.5*(sigma2_s + sigma2_l) ;
fig_number = 17;

% !!! parameters for "production" figure for paper:
Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^5 ; obs_are = 'EquallySpaced' ;
sigma2_s = 0.5 ; sigma2_l = 1.0 ; r2 = 1.0 ;
fig_number = 19;

% Nx = 25 ; Ny = Nx ; Ne = 10 ; Nreps = 10^5 ; obs_are = 'EquallySpaced' ;
% sigma2_s = 0.5 ; sigma2_l = 1.0 ; r2 = 1.0 ;
% fig_number = 21;


% !!! test ... case where l-s analysis has no skill w/o localization?
Nx = 25 ; Ny = Nx ; Ne = 4 ; Nreps = 10^5 ; obs_are = 'EquallySpaced' ;
sigma2_s = 1.0 ; sigma2_l = 0.00 ; r2 = 1.0 ;
fig_number = 99;

paramString = [ ...
        'Nx = ' num2str(Nx), ', Ny = ', num2str(Ny) ...
        ', Ne = ', num2str(Ne) ...
        ', sigma2_s = ' num2str(sigma2_s) ...
        ', sigma2_l = ' num2str(sigma2_l) ...
        ', r2 = ' num2str(r2) ] ;
disp( paramString )

localizationList = [2 4 8 12 16 24 36 50 70 90 110];
localizationList = [2 4 8 12 24 48 200 10^4 10^5];
localizationList = [2 3 4 5 6 8 10 12 24 48];
localizationList = [2 3 4 5 6 8 10 12];
localizationList = [2 4 8 12 20 30 50];

localizationList = [1 2 3 4 5 6  8  12  16  20 30 40]; % for production run
localizationList = [1 2 4  8 16 32 64 128]; % for !!! test



P = sigma2_s * eye(Nx) + sigma2_l * ones(Nx) ;
[E,DD] = eig(P) ; sqrtP = E * sqrt(DD) * E' ; 
d = diag(DD); invsqrtP = E * diag( sqrt(1./d) ) * E' ;

H = setupH(Nx,Ny,obs_are) ; %%%%%%%%%% ... or replace with H == I

Pa = P - (H*P)' * inv( H*P*H' + r2*eye(Ny) ) * (H*P); 
mseExact = mean(diag(Pa));
  % exact expected MSE from KF
mseExactAv = mean(Pa(:)) ; 
  % exact expected MSE for elementwise average
Pa_dev = (eye(Nx) - ones(Nx)/Nx) * Pa * (eye(Nx) - ones(Nx)/Nx) ; 
mseExactDev = mean(diag(Pa_dev)) ; 
  % exact expected MSE for deviations from average

Htilde = (1/sqrt(r2)) * H * sqrtP ; % H = I                                    
[U,Lambda,V] = svd(Htilde) ;

corFunHandle = @correlation_GC ; % to do periodic, only change needed is to
                                 % replace with @correlation_GCperiodic
                                 % (I think!)

mseFnLoc = zeros(length(localizationList),3) ;
mseO     = zeros(length(localizationList),3) ;
mseTransfFnLoc= zeros(length(localizationList),3) ;
mseFnLoc_test = zeros(length(localizationList),3) ;
for iloc = 1:length(localizationList)
    
localizationScale = localizationList(iloc) ; 
disp(['localizationScale = ' num2str(localizationScale)])
C = makeCorrelationMatrix( 0:Nx-1, localizationScale, 0., corFunHandle ) ;

mse_o = zeros(Nreps,3) ; mse_loc = zeros(Nreps,3) ; mse_transf = zeros(Nreps,3) ;
mse_loc_test = zeros(Nreps,3) ;
for ii = 1:Nreps
    x = sqrtP * randn(Nx,1) ; y = H*x + sqrt(r2) * randn(Ny,1);
    X = sqrtP * randn(Nx,Ne); X = (X - repmat( mean(X,2), [1 Ne] ))/sqrt(Ne);
    %X = sqrtP * randn(Nx,Ne); X = X - repmat( mean(X,2), [1 Ne] );
    xf = zeros(Nx,1) ;
    
    % optimal transf, then update
    Xtil = V' * invsqrtP * X ; % apply optimal transf to ensemble and obs
    dytil = U' * (1/sqrt(r2)) * ( y - H * xf ) ;
    Ptil = diag( diag( Xtil * Xtil' ) ) ; % "complete" localization (keep only diagonal)
    Dtil = Lambda * Ptil * Lambda' + eye(Ny) ;
    xtila = V' * invsqrtP * xf + (Ptil * Lambda') * ( Dtil \ dytil ) ;
                % update in transf variables, with localized sample cov
    xa_opt = sqrtP * V * xtila ;  % update in original variables
    
    % std localized update
    %%%Ploc = C .* P ; Dloc = Ploc + r2*eye(Ny) ;  % Test: real P, localized
    Ploc = C .* ( X * X' ) ; Dloc = H * Ploc * H' + r2*eye(Ny) ;
    %Ploc = C .* (1/Ne) * ( X * X' ) ; Dloc = H * Ploc * H' + r2*eye(Ny) ;
    xa_loc = Ploc * H' * ( Dloc \ (y - H * xf) ) ;
    K = Ploc * H * inv(Dloc) ; ImKH = eye(Nx) - K * H ;
    Pa_loc = ImKH * P * ImKH' + r2 * (K * K') ; % cov of general (possibly suboptimal) update, for prior P
    
    % transform to elementwise mean and deviations; assume ind; 
    % update bar as scalar (no loc) and devs with loc
    ybar = mean(y)  ; ydev = y - ybar ; 
    xfbar = mean(xf)  ; Xbar = mean(X,1); 
    xfdev = xf - xfbar; Xdev = X - repmat(Xbar,[Nx 1]); 
                        Pdev = C .* (Xdev * Xdev') ; Ddev = H * Pdev * H' + r2*eye(Ny) ;
    xabar = ( var(Xbar) / ( var(Xbar) + r2/Nx ) ) * ( ybar - xfbar ) ;
    xadev = Pdev * H' * ( Ddev \ (ydev - H * xfdev ) ) ;
    
    mse_transf(ii,2) = ( mean(x) - xabar )^2 ;
    mse_transf(ii,3) = mean( ( (x-mean(x)) - xadev ).^2 ) ;
    mse_transf(ii,1) = mean( ( xabar + xadev - x ).^2 ) ;
    
    mse_o(ii,1) = mean( (x - xa_opt).^2 ) ;
    mse_o(ii,2) = (mean(x) - mean(xa_opt)).^2  ;
    mse_o(ii,3) = mean( (x - xa_opt -mean(x-xa_opt)).^2 ) ;
 
    mse_loc(ii,1) = mean( (x - xa_loc).^2 ) ; 
    mse_loc(ii,2) = (mean(x) - mean(xa_loc)).^2  ;
    mse_loc(ii,3) = mean( (x - xa_loc -mean(x-xa_loc)).^2 ) ;
    
    mse_loc_test(ii,1) = mean(diag(Pa_loc)) ;

end

mseFnLoc(iloc,:) = mean( mse_loc, 1 ) ;
mseO(iloc,:) = mean( mse_o, 1 ) ;
mseTransfFnLoc(iloc,:) = mean( mse_transf, 1 ) ;
mseFnLoc_test(iloc,1) = mean( mse_loc_test(:,1) ) ;
end

mean_mseO = mean(mseO,2) ;

% Exact MSE from KF ... to compare with figs, divide by sigma2_s + sigma2_l
Pa = P - (H*P)' * inv( H*P*H' + r2*eye(Ny) ) * (H*P); mseExact = mean(diag(Pa))  
    % posterior covariance, exact KF
mean(Pa(:))   % MSE for elementwise average, exact KF
Pa_dev = (eye(Nx) - ones(Nx)/Nx) * Pa * (eye(Nx) - ones(Nx)/Nx) ; mean(diag(Pa_dev))
    % posterior cov for deviations, and MSE, exact KF


tmp_2scalesPlot