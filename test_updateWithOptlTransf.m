% Demonstrate localization after optl transf improves on "std" localization
%
% 1D example (can choose cov for state process).  Compare MSE and posterior
% variance from these updates: KF, no-loc EnKF, "std" loc EnKF, EnKF after 
% optl transf then retaining only diagonals.
%
% From test_LambdaRMSE.m ; 8 Dec 2018
% (which came from test_LinearTransformInhomogeneous.m ; 20 April 2018)
% 
% Demonstrate benefit of localizaiton after optimal LT:
%   -- specify P, H, R; draw ens from P; compute updates (see above)
%   -- loop over many realizations, and accumulate stats


%-- Initialize random number generator with specified seed -------------
rng(state_rng) ;  % ... to get same xt, y, etc w/ multiple calls to this script

%-- Some derived variables ----------------------
iloc = 1:Nx ;               % locations on 1D grid
xloc = (xDomain/Nx) * (iloc-1) ; % spatial locations

%-- Set up obs operator ----------------------
H = setupH(Nx,Ny,obs_are) ;

%-- set up covariance P --------------------
switch covIs
    case 'single'
        P = makeCorrelationMatrix( xloc, correlationScale, nu, corFunHandle ); 
    
% "multiscale versions" ... makes loc EnKF perform more poorly??
    case 'multiI'
% version I: sum of three isotropic, homogenous with different length scales 
        P = ( P + ...
   makeCorrelationMatrix( xloc, 0.2*correlationScale, nu, corFunHandle ) + ...
    makeCorrelationMatrix( xloc, 5*correlationScale, nu, corFunHandle ) ) / 3;
    case 'multiII'
% version II: small scale in left half of domain, large scale in the other
        L = 10 ; % half-width of transition between "small scale" and "large scale" 
         % parts of domain 
        a = 0.9; % proportion of variance in dominant scale in each part of domain
        variance_s = linTransition(1:Nx, Nx,L,a,1-a) ;
        variance_l = 1 - variance_s ; 
        Pl = makeCorrelationMatrix( xloc, 5*correlationScale, nu, corFunHandle ) ;
        P =   diag(variance_s) * P * diag(variance_s) ...
            + diag(variance_l) * Pl * diag(variance_l) ;
    case 'simple'
% simplest two-scale example, dating back to 2014!
        sigma2_s = 0.5 ; sigma2_l = 0.5 ;
        P = sigma2_s * eye(Nx) + sigma2_l * ones(Nx) ;
end

%-- set up localization matrix C --------------------
C = makeCorrelationMatrix( xloc, localizationScale, 0., locFunHandle) ;

%... to test (below) invsqrtC for transform (rather than invsqrtP), save these:
%[Ec,Dc]=eig(C);
%sqrtC = Ec * sqrt(Dc) * Ec' ;
%invsqrtC = Ec * inv(sqrt(Dc)) * Ec' ;

if 0 == 1; figure(1); clf; plot( xloc, P(:,Nx/2), 'kx' )
hold on; plot( xloc, P(:,Nx/4), 'kx' )
end

%-- set up obs covariance D and (correct, optimal) Pa --------------------
D = H * P * H' + r2 * eye(Ny) ;  % D = H P H^T + R ; useful in update
invD = inv(D) ;  % explicit inv likely not super accurate in big problems
                 % ... don't use in mean update below
Pa = P - P * H' * invD * H * P ;
mse_KF_exact = mean(diag(Pa)) ;

%-- set up eig(P), sqrtP, and Htilde and its SVD (which gives exact optl transform) --------------------
[E,DD] = eig(P) ; d = diag(DD) ;
                  d = (d > 0) .* d ; DD = diag(d) ;
                  % With Gaussian spatial cov, can have evals < 0
                  % (rounding errors)
sqrtP = E * sqrt(DD) * E' ;
invsqrtP = E * diag( sqrt(1./d) ) * E' ; % Note this won't work if some evals = 0!
%jnk = sqrtP * invsqrtP ; jnk(1:10,1:10)

Htilde = (1/sqrt(r2)) * H * sqrtP ;                                          
[U,Lambda,V] = svd(Htilde) ;

% figure; plot(diag(Lambda)); pause
tmp = diag(Lambda);
if Ny >= 20
    disp( 'lambdas :    1     2     3    4    5    10      20 ' )
    disp( [ tmp(1:5)' tmp(10) tmp(20) ] )
elseif Ny >= 10
    disp( 'lambdas :    1     2     3    4    5    10 ' )
    disp( [ tmp(1:5)' tmp(10) ] )
elseif Ny >= 4
    disp( 'lambdas :    1     2     3    4   ' )
    disp( tmp(1:4)' )
end
% disp( [ tmp(1:5)' tmp(10) tmp(20) ] - 10 )

%-- Do many realizations of updates ---------------------------
mse_KF = zeros(Nreps,1); mse_EnKF = zeros(Nreps,1); mse_loc = zeros(Nreps,1); mse_opt = zeros(Nreps,1) ; mse_OptLoc = zeros(Nreps,1) ;
vara_KF = zeros(Nreps,1); vara_EnKF = zeros(Nreps,1); var_loc = zeros(Nreps,1); vara_opt = zeros(Nreps,1) ;
for iloop = 1:Nreps 
    %-- generate xt, y, ensemble
    xt = sqrtP * randn(Nx,1) ; 
    y = H * xt + sqrt(r2) * randn(Ny,1) ;
    X = sqrtP * randn(Nx,Ne) ; X = (1/sqrt(Ne)) * (X - repmat( mean(X,2), [1 Ne]) );
    
    xf = zeros(Nx,1) ; % all methods will get same prior mean (exact!)
    % 7 April 2020: Try with ens mean as prior mean
    xf = mean(X,2) ;
    
    Phat = X * X' ; Dhat = H * Phat * H' + r2*eye(Ny) ;  % useful in update
    Phatloc = C.*Phat ; Dhatloc = H * Phatloc * H' + r2*eye(Ny) ;  % useful in update
    
    % updates for KF, (no loc) EnKF, localized EnKF
    xa_KF    = xf + ( P * H' ) * ( D \ (y - H * xf ) ) ;
    xa_EnKF  = xf + ( Phat * H' ) * ( Dhat \ (y - H * xf ) ) ;
    xa_loc   = xf + ( Phatloc * H' ) * ( Dhatloc \ (y - H * xf ) ) ;
    
    % optimal transf, then update
    Xtil = V' * invsqrtP * X ; % apply optimal transf to ensemble and obs
    dytil = U' * (1/sqrt(r2)) * ( y - H * xf ) ;
    Ptil = diag( diag( Xtil * Xtil' ) ) ; % "complete" localization (keep only diagonal)
%    Dtil = Lambda * Ptil * Lambda' + eye(Ny) ;
%    xtila = V' * invsqrtP * xf + (Ptil * Lambda') * ( Dtil \ dytil ) ;
    % Lambdatrunc = (Lambda > 2/sqrt(Ne)) .* Lambda ; % remove small lambdas
    Lambdatrunc = (Lambda > 0) .* Lambda ; % remove small lambdas
    Dtil = Lambdatrunc * Ptil * Lambdatrunc' + eye(Ny) ;
    xtila = V' * invsqrtP * xf + (Ptil * Lambdatrunc') * ( Dtil \ dytil ) ;
                % update in transf variables, with localized sample cov
    xa_opt = sqrtP * V * xtila ;  % update in original variables
    
    PhatOptLoc = sqrtP * V * diag( diag( Xtil * Xtil' ) ) * V' * sqrtP' ; % complete localization, then transform back
                % checked: same as updating in transf vars (above)
    DhatOptLoc = H * PhatOptLoc * H' + r2*eye(Ny) ;
    xa_opt1 = xf + ( PhatOptLoc * H' ) * ( DhatOptLoc \ (y - H * xf) ) ;
     
    % transf on state only using evecs of P, localize, then update in original variables
    XtilE = invsqrtP * X ; % this appears to be best approach (optimal?!?) 
    PhatOptLocE = sqrtP * diag( diag( XtilE * XtilE' ) ) * sqrtP' ; % complete localization, then transform back
    %XtilE = invsqrtC * X ; % try this ... still good? NO
    %PhatOptLocE = sqrtC * diag( diag( XtilE * XtilE' ) ) * sqrtC' ; % complete localization, then transform back
    DhatOptLocE = H * PhatOptLocE * H' + r2*eye(Ny) ;
    xa_OptLoc = xf + ( PhatOptLocE * H' ) * ( DhatOptLocE \ (y - H * xf) ) ;
    
    XtilE1 =  E' * X ;  % transform ensemble into space of evecs of P
    PhatOptLocE1 = E * diag( diag( XtilE1 * XtilE1' ) ) * E' ; % complete localization, then transform back
    %XtilE1 = diag( sqrt(1./d) ) * E' * X ;  % better to whiten? (vs above)
    % ... NO; mathly equiv, since diagonal scaling cancels in PhatOptLocE1 
    %PhatOptLocE1 = E * sqrt(DD) * diag( diag( XtilE1 * XtilE1' ) ) * sqrt(DD) * E' ; % complete localization, then transform back
    DhatOptLocE1 = H * PhatOptLocE1 * H' + r2*eye(Ny) ;
    xa_OptLoc1 = xf + ( PhatOptLocE1 * H' ) * ( DhatOptLocE1 \ (y - H * xf) ) ;
                                                  
    % update of covariances  
    Pa_EnKF = X * ( eye(Ne)  - (H*X)' * inv(Dhat) * (H*X) ) * X' ;
    % old: Pa_loc =  Phatloc - Phatloc * H' * inv(Dhatloc) * H * Phatloc ;
        % this is Pa from KF, in the case that P = Phatloc. It's NOT Pa
        % predicted by localized EnKF ... which is below
    Khat_loc = Phatloc * H' * inv(Dhatloc) ; ImKloc = eye(Nx) - Khat_loc * H;
    Pa_loc =  ImKloc * Phat * ImKloc' + Khat_loc * (r2*eye(Ny)) * Khat_loc' ;

    % Pa_opt = sqrtP * V * ( eye(Nx)  - Lambda' * inv(Dtil) * Lambda ) * V' * sqrtP ;
    % Pa_opt = sqrtP * V * ( Ptil  - Lambda' * inv(Dtil) * Lambda ) * V' * sqrtP ;
    % old: Pa_opt = sqrtP * V * ( Ptil  - Ptil*Lambda' * inv(Dtil) * Lambda*Ptil ) * V' * sqrtP ;
        % like above, this is Pa in case that P = Ptil, which is NOT the
        % same as Pa predicted by localized EnKF in opt coords, below:    
    Khat_loc = Ptil * Lambdatrunc' * inv(Dtil) ; ImKloc = eye(Nx) - Khat_loc * Lambdatrunc;
    Pa_opt =  sqrtP * V * ( ...
        ImKloc *  (Xtil * Xtil') * ImKloc' + Khat_loc * Khat_loc' ...
                          ) * V' * sqrtP ;

    % compute mse in opt coords, if desired, rather than origl state space
    % TO DO: Also transform vara's !!!
    if 0 == 1
        Tmp = V' * invsqrtP;
        xt = Tmp*xt         ; xa_KF  = Tmp*xa_KF  ; xa_EnKF   = Tmp*xa_EnKF ; 
        xa_loc = Tmp*xa_loc ; xa_opt = Tmp*xa_opt ; xa_OptLoc = Tmp*xa_OptLoc ; 
    end
    
    % posterior mse and variance
    mse_KF(iloop)  = mean((xt - xa_KF).^2) ;
    mse_EnKF(iloop)= mean((xt - xa_EnKF).^2) ;
    mse_loc(iloop) = mean((xt - xa_loc).^2) ;
    mse_opt(iloop) = mean((xt - xa_opt).^2) ;
    mse_opt1(iloop) = mean((xt - xa_opt1).^2) ;
    mse_OptLoc(iloop) = mean((xt - xa_OptLoc).^2) ;
    mse_OptLoc1(iloop) = mean((xt - xa_OptLoc1).^2) ;
    
    vara_KF(iloop)  = mse_KF_exact ;
    vara_EnKF(iloop)= mean( diag(Pa_EnKF) ) ;
    vara_loc(iloop) = mean( diag(Pa_loc ) ) ;
    vara_opt(iloop) = mean( diag(Pa_opt ) ) ;

end  

