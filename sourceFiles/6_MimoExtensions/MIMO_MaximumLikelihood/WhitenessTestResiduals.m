function [Auto_Corr, Lags, Conf_Bound, Fraction] = WhitenessTestResiduals(G, varG, Gest, CL, dof, PlotFig);
%
%		Whiteness test on the difference between the measured G (nonparametric estimate) and modeled Gest 
%       (parametric estimate) frequency response matrix, taking into account the possible correlation over the 
%       frequency (CL) and the variability (dof) of the noise covariance estimate varG. 
%
%       function [Auto_Corr, Lags, Conf_Bound] = WhitenessTestResiduals(G, varG, Gest, CL, dof, PlotFig);
%
%       References:
%
%                   Pintelon, R., G. Vandersteen, J. Schoukens, and Y. Rolain (2011). Improved (non-)parametric identification of dynamic 
%                   systems excited by periodic signals - The multivariate case, Mechanical Systems and Signal Processing, vol. 25, no. 8, 
%                   pp. 2892-2922.
%
%                   Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%                   IEEE Press-Wiley, Piscataway (USA). 
%
%
%   Output
%
%       Auto_Corr   =   autocorrelation of each entry of the residuals (G-Gest)/std(G) 
%                           size  ny x nu x 2*floor(F/(CL+1))+1 
%
%       Lags        =   lags (in multiples of CL+1) over which the auto-correlation is calculated  
%
%       Conf_Bound  =   50% (Conf_Bound(1,:)) and 95% (Conf_Bound(2,:)) confidence bound auto-correlation
%                           size 2 x 2*floor(F/(CL+1))+1 
%
%       Fractions   =   fraction in % outside the 50% and 95% confidence bounds
%                           size ny x nu x 2
%
%       
%   Input
%
%       G           =   measured frequency response matrix (FRM)
%                           size ny x nu x F
%
%       varG       =    either: variance of each entry of the measured FRM G 
%                           size ny x nu x F
%                       or: covariance of vec(G)
%                           size ny*nu x ny*nu x F
%
%       Gest        =   estimated plant transfer function model 
%                           size ny x nu x F
%
%       CL          =   ± CL is the correlation length over the frequency of the measured FRM G 
%                       (optional, default = 0 = no correlation)
%
%       dof         =   degrees of freedom of the variance estimate varG, optional parameter
%                       (optional, default = inf)
%
%       PlotFig     =   plot auto-correlation if 1, else no figures
%                       (optional, default = 1)
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 19 November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 10 October 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    if isempty(CL)
        CL = 0;
    end % if   
catch
	CL = 0;    
end % try

[ny, nu, F] = size(G);
Res = zeros(ny, nu, F);
% F1 = number of independent FRM measurements 
% if no correlation (CL = 0), then F1 = F
F1 = floor((F-1)/(CL+1))+1;       
Auto_Corr = zeros(ny, nu, 2*F1-1);
Fraction = zeros(ny, nu, 2);
Lags = [-F1+1:1:F1-1];

try
    if isempty(dof)
        dof = inf;
    end % if   
catch
	dof = inf;    
end % try

try
    if isempty(PlotFig)
        PlotFig = 1;
    end % if   
catch
	PlotFig = 1;    
end % try

% check whether variance of each entry FRM or covariance vec(G) is provided
rows = size(varG, 1);
columns = size(varG, 2);

if (rows-ny+columns-nu) > 0
    % keep the diagonal of CvecG only
    dummy = varG;
    varG = zeros(ny,nu,F);
    for kk = 1:F
        varG(:,:,kk) = reshape(diag(dummy(:,:,kk)), [ny, nu]);
    end % kk
    clear dummy
end % if cov(vec(G)) is provided


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of the auto-correlation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalised residuals
for kk = 1:F    
    Res(:,:,kk) = (G(:,:,kk) - Gest(:,:,kk))./varG(:,:,kk).^0.5;    
end % kk, frequencies

% auto-correlation
for ii = 1:nu
    for jj = 1:ny
        Auto_Corr(jj,ii,:) = xcorr(squeeze(Res(jj,ii,1:CL+1:end)),'unbiased');
    end % jj, outputs
end % ii, inputs

% scaling auto-correlation to account for the degrees of freedom of the
% variance estimate varG 
switch isinf(dof)
    case 0
        scale0 = (dof-1)/dof;
        scale = (dof-2/3)/(dof+1/12);
        % factor 3 for noise models from LocalPoly and FastLocalPoly
        % for RobustLocalPoly => 3 is an upperbound
        Conf_scale0 = scale0*sqrt(3*dof^3/(dof-1)^2/(dof-2));         
        Conf_scale = scale*dof/(dof-1);
    case 1
        scale0 = 1;
        scale = 1;
        Conf_scale0 = 1;
        Conf_scale = 1;
end % switch
TheScale = scale*ones(size(Lags));
TheScale(F1) = scale0;
Auto_Corr = Auto_Corr .* repmat(reshape(TheScale, [1,1,2*F1-1]), [ny, nu, 1]);

% std auto-correlation
p50 = sqrt(-log(1-0.5));
p95 = sqrt(-log(1-0.95));
ConfScale = Conf_scale*ones(size(Lags));
ConfScale(F1) = Conf_scale0;
Conf_Bound = repmat(ConfScale./(F1-abs(Lags)).^0.5, [2,1]);
Conf_Bound(1,:) = p50*Conf_Bound(1,:);
Conf_Bound(2,:) = p95*Conf_Bound(2,:);

% fractions outside p% confidence bounds
Select = [1:1:2*F1-1];
Select(F1) = [];
for ii = 1:nu
    for jj = 1:ny
        Fraction(jj,ii,1) = length(find(squeeze(abs(Auto_Corr(jj,ii,Select))) - Conf_Bound(1,Select).' > 0))/(2*F1-2)*100;
        Fraction(jj,ii,2) = length(find(squeeze(abs(Auto_Corr(jj,ii,Select))) - Conf_Bound(2,Select).' > 0))/(2*F1-2)*100;
    end % jj, outputs
end % ii, inputs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the auto-correlation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PlotFig == 1
    figure
    mm = 0;
    for jj = 1:ny
        for ii = 1:nu
            mm = mm+1;
            subplot(ny,nu,mm);
            plot(Lags, squeeze(abs(Auto_Corr(jj,ii,:))), '*', Lags, Conf_Bound(1,:), 'k', ...
                 Lags, Conf_Bound(2,:), 'k--') 
            xlabel('{\itk}')
            Rlabel = ['[',num2str(jj),',',num2str(ii),']'];
            ylabel(['{\itR}_{',Rlabel,'}({\itk})'])
        end % jj, outputs
    end % ii, inputs
    zoom on
    shg
end % if auto-correlation should be plotted

