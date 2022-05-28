% THIS VERSION: finds best fixed localization scale as fn of Ne;
%   see use of localizationList below. 
% 11 March 2022: Seems to be only a start; definitely doesn't work.
%               ... modifying to follow wrapper_updateAsFnAnything
%
% wrapper for test_updateWithOptTransf.m; began from wrapper_updateAsFnNe.m
%
% Various input parameters are controlled here (rather than buried in 
% test_updateWithOptlTransf.m).  Probably should pull more of setup out of
% that script too.  List of parameters that are here (10 May 2019):
%   xDomain, Nx, Ny, obs_are, r2, correlationScale, corFunHandle, nu,
%   localization Scale, localizationAsFnNe, locFunHandle, Ne, Nreps,
%   paramName, paramList
%
% Does loop over a chosen parameter; saves/plots 

rng('shuffle'); state_rng = rng; 
            % Get new seed/state for random number generator; save for use
            % in each realization in loop

xDomain = 2*pi ;  % width of spatial domain [why not == 1??]
Nx = 200 ;  % number of spatial points
Ny = Nx / 4 ; 
            % number of obs
            % [Require Nx/Ny = integer for 'EquallySpaced' ???]
obs_are = 'EquallySpaced' ; 
            % kind of obs network: 
            %   EquallySpaced == pt obs, equally spaced, iid errors
            %   Random        == pt obs, randomly located, iid errors
            %   Averages      == local average of state (to be implemented)
r2 = 10^0 ; % obs-error variance (relative to prior variance, == 1 )           
correlationScale = 0.05 * xDomain / sqrt(2);  % 22 Nov 2019: factor of 1/sqrt(2) with new version of correlation_matern.m
            % prior length scale in units of domain size
covIs = 'single' ;
            % single-scale or (various) multiscale covariances
corFunHandle = @correlation_matern ; 
            % correlation fn for prior covariance; Matern and sqd expl work
% >> could also try correlation specifed by spectral power law
%  ... see comments in test_LinearTransformInhomogeneous.m for subtleties
% >> other complexification of prior cov: superpose multiple correlation
% fns with difft scales; inhomogeneous variance and/or correlation str.
nu = 3/2 ;  % Matern "smoothness" parameter; sample paths have nu - 1 deri
localizationScale = 2*correlationScale ;  
            % localization scale in units of domain size [check dependence!!!]
            % ... if switched to Matern, need to modify call to
            %     makeCorrelationMatrix in test_updateWithOptlTransf.m
locFunHandle = @correlation_GC ;   
            % correlation fn for localization; Matern, sqd expl, GC
Ne = 4 ;    % ensemble size
Nreps = 10^4 ; % used for results in mmse_6April_Ne1024.mat
           % rmse etc calculated over this many realizations
Nreps = 10^2 ; 

% Parameter to vary
paramName = 'ensemble size' ;
paramList = 2.^(2:10) ; % used for results in mmse_6April_Ne1024.mat
paramList = 2.^(2:8) ; 

localizationList = [16 8 4 2 1 0.5 0.37 0.25 0.19 0.125] * correlationScale ;
localizationList = [16 10 8 6 5 4.5 4 3.5 3 2] * correlationScale ;

% !!! TO CHANGE PARAMETER THAT VARIES, change assignment at beginning of
% "iwrapper" loop below:  xxxx = paramList(iwrapper)
if 1 == 1
% arrays to store results on "grid"
mmse_EnKF   = zeros( length(paramList), length(localizationList) ) ;
mmse_loc    = zeros( length(paramList), length(localizationList) ) ;
mmse_optloc = zeros( length(paramList), length(localizationList) ) ;
mvara_EnKF   = zeros( length(paramList), length(localizationList) ) ;
mvara_loc    = zeros( length(paramList), length(localizationList) ) ;
mvara_optloc = zeros( length(paramList), length(localizationList) ) ;
 
for iwrapper = 1:length(paramList)
    %correlationScale = paramList(iwrapper)
    %localizationScale = 2*correlationScale ;
    Ne = paramList(iwrapper) 
    for jlocalization = 1:length(localizationList) 
        localizationScale = localizationList(jlocalization) 
        test_updateWithOptlTransf
        mmse_EnKF  (iwrapper,jlocalization) = mean(mse_EnKF) ;
        mmse_loc   (iwrapper,jlocalization) = mean(mse_loc) ;
        mmse_optloc(iwrapper,jlocalization) = mean(mse_opt) ;
        varmse(iwrapper,jlocalization) = var(mse_loc) ;
        mvara_EnKF  (iwrapper,jlocalization) = mean(vara_EnKF) ;
        mvara_loc   (iwrapper,jlocalization) = mean(vara_loc) ;
        mvara_optloc(iwrapper,jlocalization) = mean(vara_opt) ;
    end
%    mse(:,iwrapper) = [mean(mse_KF), mean(mse_EnKF), mean(mse_loc), mean(mse_opt) mean(mse_OptLoc) ] ;
%    vara(:,iwrapper) = [mean(vara_KF), mean(vara_EnKF), mean(vara_loc), mean(vara_opt) ] ;
%    vara_exact_KF(iwrapper) = mse_KF_exact ;
%    disp( [mse_KF_exact mean(vara_KF) mean(vara_KF)/mse_KF_exact])
end
end

return
