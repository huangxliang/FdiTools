function [CZ, Z, freq, G, CvecG, dof, CL] = RobustLocalPolyAnal(data, method);
%
%       Using P consecutive periods and M independent realisations of nu independent MIMO experiments 
%       on an nu input and ny output system, the following quantities are calculated:
%
%           (i)     The generalized sample means, sample (co-)variances of the nu input/output DFT spectra 
%                   with optimal noise and/or system transient suppression,
%
%           (ii)	The nu x ny frequency response matrix (FRM) and its covariance. 
%
%           (iii)	The covariances are split in the contribution of the noise and the stochastic nonlinear distortions. 
%       
%       Notes:
%
%           (i)     At least P = 2 periods should be measured for calculating the noise covariances. 
%
%           (ii)    At least M = 2 realisations should be measured for calculating the total variance (noise + 
%                   stochastic nonlinear distortions).
%
%           (iii)   At least P = 2 periods and M = 2 realisations are needed to estimate the stochastic nonlinear distortions. 
%
%           (iv)    If the generator noise is dominant, then CvecG.NL >> CvecG.NLcheck and, the averaging over the realisations 
%                   should be done on the FRM's per realisation (no bias is introduced). This is done automatically if no 
%                   reference signal is provided (data.r is empty). The generator noise is typically dominant over the multisine
%                   realisations for lowly damped mechanical systems excited by a shaker.
%
%           (v)     A local polynomial approximation of the noise (and system) transient term is made, while the FRM is NOT 
%                   approximated.
%
%       References:
%
%                   Pintelon, R., K. Barbé,  G. Vandersteen, and J. Schoukens (2011). Improved (non-)parametric identification 
%                   of dynamic systems excited by periodic signals, Mechanical Systems and Signal Processing, 
%                   vol. 25, no. 7, pp. 2683-2704. 
%
%                   Pintelon, R., G. Vandersteen, J. Schoukens, and Y. Rolain (2011). Improved (non-)parametric identification 
%                   of dynamic systems excited by periodic signals - The multivariate case, Mechanical Systems 
%                   and Signal Processing, vol. 25, no. 8, pp. 2892-2922.
%
%                   Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%                   IEEE Press-Wiley, Piscataway (USA).
%
%       function [CZ, Z, freq, G, CvecG, dof, CL] = RobustLocalPolyAnal(data, method);
%
%
%   Output
%
%       CZ      =   sample covariance matrices of the sample mean of the Fourier coefficients Z for the nu independent MIMO 
%                   experiments, where Z = [Y; U] contains the ouput and input spectra stacked on top of each other; 
%                   size (ny+nu) x (ny+nu) x nu x F 
%                   struct('n', [], 'NL', [], 'S', [])
%                       CZ.n            =   noise covariance matrix of Z as complex numbers 
%                                               CZ.n(:, :, ii) = Cov(NZ) /(M*P)
%                                               CZ.n(1:ny, 1:ny, ii) = Cov(NY) /(M*P)
%                                               CZ.n(ny+1:end, ny+1:end, ii) = Cov(NU)/(M*P)
%                                               CZ.n(1:ny, ny+1:end, ii) = Cov(NY, NU)/(M*P) 
%                                               CZ.n(ny+1:end, 1:ny, ii) = Cov(NU, NY)/(M*P) 
%                                           with NX the noise on X w.r.t. one period and one MIMO experiment 
%
%                       CZ.NL           =   total (noise + NL distortion) covariance matrix of Z 
%                                               CZ.NL(:, :, ii) = Cov(NZ) /(M*P) + Cov(ZS) / M 
%                                               CZ.NL(1:ny, 1:ny, ii) = Cov(NY) /(M*P) + Cov(YS) / M 
%                                               CZ.NL(ny+1:end, ny+1:end, ii) = Cov(NU) /(M*P) + Cov(US) / M 
%                                               CZ.NL(1:ny, ny+1:end, ii) = Cov(NY, NU)/(M*P) + Cov(YS, US) / M  
%                                               CZ.NL(ny+1:end, 1:ny, ii) = Cov(NU, NY)/(M*P) + Cov(US, YS) / M  
%                                           with NX the noise on X w.r.t one period and one MIMO experiment, and 
%                                           with XS the stochastic NL distortions w.r.t. one realisation and one MIMO experiment 
%
%                       CZ.S            =   covariance matrix of the stochastic nonlinear distortions in Z w.r.t. one realisation,
%                                           and one MIMO experiment
%                                               CZ.S(:, :, ii) = Cov(ZS) 
%                                               CZ.S(1:ny, 1:ny, ii) = Cov(YS) 
%                                               CZ.S(ny+1:end, ny+1:end, ii) = Cov(US) 
%                                               CZ.S(1:ny, ny+1:end, ii) = Cov(YS, US)  
%                                               CZ.S(ny+1:end, 1:ny, ii) = Cov(US, YS)  
%                                           with XS the stochastic NL distortions w.r.t. one realisation and one MIMO experiment
%
%                   Notes: 
%                        - CZ.NL, and CZ.S are calculated only if y(t), u(t), and r(t) are available 
%                        - if the actuator is linear and there is no nonlinear interaction between the generator and the
%                          NL system, then US = 0 
%
%       Z       =   sample mean Fourier coefficients with noise and/or system transient removed for all M independent full MIMO experiments 
%                   (= nu MIMO experiments), where Z stacks the output and input spectra stacked on top of each other: Z = [Y; U]; 
%                   size (ny+nu) x nu x F 
%                   If the reference signal r(t) is available then
%                       Z is the mean value of the projected output-input spectra over all periods and realisations
%                   else
%                   	Z is empty except if M = 1 then Z is the mean value over all periods
%
%       freq    =   frequency of the excited harmonics; size 1 x F 
%
%	    G		=	estimated frequency response matrix; size ny x nu x F
%
%       CvecG   =   covariance matrices of vec(G); size (ny*nu) x (ny*nu) x F 
%                   struct('n', [], 'NL', [], 'NLcheck', [], 'S')
%                       CvecG.n         =   noise covariance matrix of vec(G) 
%
%                       CvecG.NL        =   total (noise + nonlinear distortion) covariance matrix of vec(G) 
%
%                       CvecG.NLcheck   =   total (noise + nonlinear distortion) covariance matrix of vec(G) calculated
%                                           via the classical averaging over the realisations (no averaging over neighbouring 
%                                           excited frequencies) 
%
%                       CvecG.S         =   covariance matrix of the stochastic nonlinear distortions in vec(G) w.r.t. one
%                                           realisation and one MIMO experiment 
%
%       dof     =   actual degrees of freedom of the input-output and FRF covariance estimates = number of equivalent 
%                   independent experiments - 1 
%                   struct('n', [], 'NL', [], 'NLcheck', [], 'Gn', [], 'GNL', [], 'GNLcheck', [])
%                       dof.n           =   actual degrees of freedom of the noise covariances (to reduce the uncertainty 
%                                           of the leakage suppression, the minimum dof per realisation is ny + 5)  
%                       dof.NL          =   actual degrees of freedom of the total covariance 
%                       dof.Gn          =   actual degrees of freedom of the noise covariance of vec(G) 
%                       dof.GNL         =   actual degrees of freedom of the total covariance of vec(G) 
%                       dof.GNLcheck    =   actual degrees of freedom of the total covariance as calculated in CvecG.NLcheck 
%
%       CL      =   ± CL is the correlation length (over the frequency) of the sample mean and/or sample covariance 
%                   struct('n', [], 'NL', [])
%                       CL.n            =   correlation length Z, G, CZ.n, and CvecG.n over the EXCITED frequency 
%                                           value at frequency k is correlated with that at frequency m in [k-CL.n, k+CL.n]
%                       CL.NL           =   correlation length CZ.NL, and CvecG.NL over the EXCITED frequencies 
%                                           value at frequency k is correlated with that at frequency m in [k-CL.NL, k+CL.NL]
%                   Note: CL+1 = the frequency width (in samples) of the local polynomial estimate 
%
%
%   Input
%
%       data    =   structure containing the input/output time domain signals 
%                   struct('u', [], 'y', [], 'r', [], 'N', [], 'Ts', [], 'ExcitedHarm', []) 
%                       data.u              =   input signal; size nu x nu x M x N*P 
%                                               (nu inputs, nu experiments, M realisations, P periods, N points per period) 
%                       data.y              =   ouput signal; size ny x nu x M x N*P 
%                                               (ny outputs, nu experiments, M realisations, P periods, N points per period) 
%                       data.r              =   reference signal, e.g. signal stored in the arbitrary waveform generator (optional) 
%                                               size nu x nu x M x N
%                                               (nu inputs, nu experiments, M realisations, 1 period, N points per period) 
%                                               Note: - the phase spectrum varies randomly over the realisations, the experiments,
%                                                       and the frequencies, while the amplitude spectrum remains the same. 
%                                                     - if generator is dominant, then the input-output noise errors are 
%                                                       highly correlated, and the averaging over the realsisations should be done 
%                                                       on the FRM and not on the input-output spectra. This is often the case for
%                                                       lowly damped mechanical systems => do not provide the reference signal 
%                       data.N              =   number of time domain samples in one period 
%                       data.Ts             =   sampling period
%                       data.ExcitedHarm    =   excited harmonics (in harmonic numbers); size 1 x F (F = number of harmonics) 
%
%       method  =   structure containing the parameters of the method used (optional) 
%                   struct('order', [], 'dof', []) 
%                       method.order        =   order of the local polynomial approximation (default 2) 
%                       method.dof          =   degrees of freedom of the (co-)variance estimates = equivalent number of 
%                                               independent experiments - 1 (default ny) 
%                       method.transient    = 	determines the estimation of the transient term (optional; default 1)  
%                                                   1: transient term is estimated 
%                                                   0: no transient term is estimated 
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 29 October 2009 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 29 March 2010
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = data.N;                             % number of samples per period 
[ny, nu, M, NP] = size(data.y);         % nu inputs, ny outputs, and M realisations
P = NP/N;                               % number of periods
fs = 1/data.Ts;                         % sampling frequency
ExcitedHarm = data.ExcitedHarm(:).';    % excited harmonic numbers
freq = ExcitedHarm*fs/N;                % frequency excited harmonics
F = length(ExcitedHarm);                % number of excited frequencies 
dof_min = ny;                           % minimal number of degrees of freedom

try
    if isempty(method)
        method = struct('order', 2, 'dof', dof_min);
    end
catch
    method = struct('order', 2, 'dof', dof_min);
end % try

try
	if isempty(method.order)
		method.order = 2;
	end
catch
	method.order = 2;
end % try

try
	if isempty(method.dof)
		method.dof = dof_min;
	end
catch
	method.dof = dof_min;
end % try
if method.dof < dof_min
    method.dof = dof_min;
end % if

try
	if isempty(method.transient)
		method.transient = 1;
	end
catch
	method.transient = 1;
end % try
if method.transient ~= 1
    method.transient = 0;
end % if

try
    if isempty(data.r)
        Reference = 0;
    else
        Reference = 1;
    end % if
catch
    data.r = [];
    Reference = 0;
end % try

CZ = struct('n', [], 'NL', [], 'S', []);
G = zeros(ny, nu, F);
CvecG = struct('n', [], 'NL', [], 'NLcheck', [], 'S', []);
dof = struct('n', [], 'NL', [], 'Gn', [], 'GNL', [], 'GNLcheck', []);
CL = struct('n', [], 'NL', []);

if P > 1
    CZ.n = zeros(ny+nu, ny+nu, nu, F);
    CvecG.n = zeros(ny*nu, ny*nu, F);
end % if P > 1

if M > 1
    CvecG.NL = zeros(ny*nu, ny*nu, F);
    CvecG.NLcheck = zeros(ny*nu, ny*nu, F);
else % M = 1
    Z = zeros(ny+nu, nu, F);                    % mean value over the periods 
end % if M > 1

if (M > 1) && (P > 1)
    CvecG.S = zeros(ny*nu, ny*nu, F);
    if Reference
        CZ.NL = zeros(ny+nu, ny+nu, nu, F);
        CZ_NLcheck = zeros(ny+nu, ny+nu, nu, F);
        CZ.S = zeros(ny+nu, ny+nu, nu, F);  
        Z = zeros(ny+nu, nu, F);
    end % if reference signal available
end % if


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample means and sample noise (co-)variances over the P periods %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample means over the periods
SelectExcited = P*ExcitedHarm + 1;                              % selection of the excited harmonics
Zfft = fft([data.y; data.u], [], 4)/(N*P);                      % scaling DTF s.t. the correct Fourier coefficients are recovered
Z_all = Zfft(:, :, :, SelectExcited);                           % input-output Fourier coefficients

Z = zeros(ny+nu, nu, F);
% if only one period and one realisation are available
if (P == 1) && (M == 1)
    Z(:, :, :) = squeeze(Z_all);
    G = FRM(Z);
end % if one period, one realisation

% calculate the noise (co-)variances if more than one period is available
if P > 1
    
    % select all DFT lines (non-excited + excited) around the excited frequency band
    StartDFT = P*ExcitedHarm(1) - (P-1);                        % at least P-1 non-excited DFT lines exist before the first harmonic
    StopDFT = P*ExcitedHarm(end) + (P-1);                       % at least P-1 non-excited DFT lines exist after the last harmonic
    SelectAll = [StartDFT:1:StopDFT] + 1; 
    Zres_all = Zfft(:, :, :, SelectAll);                        % noisy input-output DFT lines
    
    % calculate for each experiment and each realisation the sample mean and sample 
    % noise (co-)variances of the Fourier coefficients via the local polynomial approach 
    data_IO = struct('Y', [], 'ExcitedHarm', ExcitedHarm);
    method_IO = struct('order', [], 'dof', [], 'period', P, 'transient', []);
    method_IO.order = method.order;
    method_IO.dof = ceil(method.dof/M);  
    if ~Reference
        method_IO.dof =ceil(method.dof/M/nu);                   % because of averaging covariance matrices over nu experiments
    end % if no reference signal available
    method_IO.transient = method.transient;
    CZ_all = zeros(ny+nu, ny+nu, nu, M, F);                     % intermediate variable
    
    for ii = 1:nu % loop over experiments
        for mm = 1:M % loop over realisations
            data_IO.Y = squeeze(Zres_all(:, ii, mm, :));
            [dummy, CZm, Zm, dofn, CLn] = PeriodLocalPolyAnal(data_IO, method_IO); 
            Z_all(:, ii, mm, :) = Zm;
            CZ_all(:, :, ii, mm, :) = CZm;
        end % mm, realisations
    end % ii, inputs

    CZ.n(:, :, :, :) = squeeze(mean(CZ_all, 4))/M;              % mean value over the realisations; noise covariance sample mean
    dof.n = M*dofn;                                             % actual degrees of freedom of the noise (co-)variance estimate
    CL.n = CLn;                                                 % correlation length Z and CZ.n over the frequency
    
    if M == 1
        % input-output sample means
        Z = zeros(ny+nu, nu, F);
        Z(:, :, :) = squeeze(Z_all);
        % FRM and inverse input matrix
        [G, SUinv] = FRM(Z);
        % noise covariance matrix FRM
        CvecG.n = Cov_FRM(G, SUinv, CZ.n);
        dof.Gn = nu*dof.n;
    end % if M = 1
    
end % P > 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample means and sample noise (co-)variances over the M realisations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the reference signal is available
if (M > 1) && Reference
    
    % DFT spectra of the reference signal for all realisations and all experiments 
    R_all = fft(data.r, [], 4)/N;                               % only one period of the reference signal
    R_all = R_all(:, :, :, ExcitedHarm + 1);                    % selection of the excited harmonics
    Tphase_all = R_all./abs(R_all)/sqrt(nu);                    % unitary matrix of the input phases
        
    % turning the phases of the input-output spectra for all realisations and all frequencies 
    if nu == 1
        Z_all = Z_all .* repmat(conj(Tphase_all), [ny+nu, 1, 1, 1]);
    else % nu > 1
        for mm = 1:M
            for kk = 1:F
                Z_all(:, :, mm, kk) = Z_all(:, :, mm, kk) * squeeze(Tphase_all(:, :, mm, kk))';
            end % kk, frequencies
        end % mm, realisations
    end % if one input
    
    % sample mean input-output spectra over the M realisations
    Zmean = mean(Z_all, 3);
    
    % input-output residuals
    Zres_all = Z_all - repmat(Zmean, [1, 1, M, 1]);
    
    % remove the singular realisation dimension
    Z(:,:,:) = squeeze(Zmean);
        
    % local polynomial estimate total covariance of the residuals over the M realisations and the nu experiments     
    data_res = struct('U', [], 'Y', [], 'freq', []);
    data_res.freq = ExcitedHarm;
    method_res = struct('order', [], 'dof', []);
    method_res.order = method.order;
    method_res.dof = ceil(method.dof/(M-1));                        % M-1 independent residuals
    method_res.transient = 0;                                       % the stochastic NL distortions are periodic (=> no leakage)
                                                                    % and the noise leakage has been removed  
    
    CZ_all = zeros(ny+nu, ny+nu, nu, M, F);                         % intermediate variable
    for ee = 1:nu
        for mm = 1:M
            data_res.Y = squeeze(Zres_all(:, ee, mm, :));
            [CC, xx1, xx2, xx3, xx4, dofNL, CLNL] = LocalPolyAnal(data_res, method_res);
            CZ_all(:, :, ee, mm, :) = CC.n;        
        end % mm, realisations
    end % ee, experiments
    clear xx1 xx2 xx3 xx4 CC

    CZ.NL(:,:,:,:) = squeeze(mean(CZ_all, 4))/(M-1);                % total sample covar sample mean. Note: (M-1) independent residuals
    dof.NL = (M-1)*dofNL;                                           % actual degrees of freedom of the total covariance estimate
    CL.NL = CLNL;                                                   % correlation length CZ.NL over the frequency
            
    % input-output stochastic NL distortions w.r.t. one realisation  
    if P > 1
        CZ.S = EstimStochNLdist(CZ.NL, CZ.n, M);
    end % if P > 1
    
    % calculation mean FRF, its total and noise covariance, and the stochastic NL distortions 
    [G, SUInv] = FRM(Z);                                            % FRM and inverse input matrix
    CvecG.NL = Cov_FRM(G, SUInv, CZ.NL);                            % total covariance matrix FRM   
    dof.GNL = nu*dof.NL;
    
    if P > 1
        CvecG.n = Cov_FRM(G, SUInv, CZ.n);                          % noise covariance matrix FRM  
        dof.Gn = nu*dof.n;
        CvecG.S = EstimStochNLdist(CvecG.NL, CvecG.n, M);
    end % P > 1
        
    % total covariances as sample covariance over the realisations at the excited harmonics
    % used for checking that the generator noise is dominant; if so then CvecG.NLcheck < CvecG.NL 
    Zres_all = permute(Zres_all, [3 1 2 4]);
    for kk = 1:F
        for ee = 1:nu
            CZ_NLcheck(:, :, ee, kk) = cov(squeeze(Zres_all(:, :, ee, kk)));
        end % ee, experiments
    end % kk, frequencies
    CZ_NLcheck = CZ_NLcheck/M;                                      % total covariance of the mean value
    dof_NLcheck = (M-1);                                            % actual degrees of freedom of the total covariance estimate
    CvecG.NLcheck = Cov_FRM(G, SUInv, CZ_NLcheck);                  % total covariance matrix FRM    
    dof.GNLcheck = nu*dof_NLcheck;                                  % actual degrees of freedom of the total covariance estimate
        
end % if M > 1 and r(t) is available


% the reference signal is not available
if (M > 1) && ~Reference

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mean FRF and noise covariance %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % FRM for all realisations
    G_all = zeros(ny, nu, M, F);
    Zm = zeros(ny+nu, nu, F);                                       % value input-output spectra of the mth realisation
    for mm = 1:M
        Zm(:,:,:) = squeeze(Z_all(:,:,mm,:)); 
        Gm = FRM(Zm);                                               % FRM and inverse input matrix
        G_all(:,:,mm,:) = Gm;
    end % mm, realisations
  
    % mean (over the realisations) estimate FRM
    G(:,:,:) = squeeze(mean(G_all, 3));
    
    % mean (over the realisations) estimate noise covariance matrix
    if P > 1
        % inverse of the mean input matrix
        SUinv = InvMeanInput(Z_all(ny+1:end,:,:,:));
        CZm = zeros(ny+nu, ny+nu, nu, F);
        CZm(:,:,:,:) = squeeze(mean(CZ_all,4))/M;                   % mean noise covariance matrix input-output spectra
        CvecG.n = Cov_FRM(G, SUinv, CZm);                           % noise covariance matrix FRM
        dof.Gn = nu*dof.n;                                          % actual degrees of freedom noise covariance estimate FRM 
                                                                    % over the M realisations and the nu experiments 
    end % if more than one period
        
    % local polynomial estimate total variance FRM
    Gres_all = G_all - repmat(reshape(G, [ny,nu,1,F]), [1,1,M,1]);  % FRM residuals; G is an ny x nu x F matrix
    data_res = struct('U', [], 'Y', [], 'freq', []);
    data_res.freq = ExcitedHarm;
    method_res = struct('order', [], 'dof', []);
    method_res.order = method.order;
    method_res.dof = ceil(method.dof/(M-1));                         % M-1 independent residuals;
    method_res.transient = 0;                                        % the stochastic NL distortions are periodic => no leakage
    CvecG_all = zeros(ny*nu, ny*nu, M, F);                           % intermediate variable
    for mm = 1:M
        data_res.Y = reshape(Gres_all(:,:,mm,:), [ny*nu, F]);
        [CC, xx1, xx2, xx3, xx4, dofNL, CLNL] = LocalPolyAnal(data_res, method_res);
        CvecG_all(:, :, mm, :) = CC.n;        
    end % mm, realisations
    clear xx1 xx2 xx3 xx4

    CvecG.NL(:,:,:) = squeeze(mean(CvecG_all, 3))/(M-1);            % total sample var sample mean. Note: (M-1) independent residuals
    dof.GNL = (M-1)*dofNL;                                          % actual degrees of freedom of the total covariance estimate
                                                                    % over the M realisations  
    CL.NL = CLNL;                                                   % correlation length CZ.NL over the frequency
    
    % estimate stochastic NL distortions on FRM w.r.t. one realisation
    if P > 1
        CvecG.S = EstimStochNLdist(CvecG.NL, CvecG.n, M);
    end % if P > 1
        
end % if M > 1 and r(t) is not available
