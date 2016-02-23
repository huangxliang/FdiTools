function [CZ, Z, freq, G, CvecG, dof, CL] = FastLocalPolyAnal(data, method);
%
%       Using P consecutive periods of 1 MIMO experiment on an nu input and ny output system, 
%       the following quantities are calculated:
%
%           (i)     The generalized sample means, and the sample covariances (noise, total, and stochastic nonlinear distortions)  
%                   of the input/output DFT spectra with optimal noise and/or system transient suppression (P > 1),
%
%           (ii)	The nu x ny frequency response matrix (FRM), its noise covariance (P > 1), its total covariance (P >= 1), 
%                   and the covariance of the stochastic nonlinear distortions (P > 1). 
%       
%       Notes:
%
%           (i)     To split the covariances in contributions of the noise and the stochastic nonlinear distortions
%                   a random phase multisine experiment is required and P > 1. In addition, for the input/output spectra, 
%                   the reference signal should be available. If not, then the input-ouput noise and stochastic nonlinear
%                   contributions are reflected to the output in CZ.NL. 
%
%           (ii)    A local polynomial approximation of the FRM and the noise (and system) transient term is made. 
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
%       function [CZ, Z, freq, G, CvecG, dof, CL] = FastLocalPolyAnal(data, method);
%
%
%   Output
%
%       CZ      =   sample covariance matrix of the sample mean of the Fourier coefficients Z, where Z = [Y; U] contains 
%                   the ouput and input spectra stacked on top of each other;  
%                   struct('n', [], 'NL', [], 'S', [], 'm_NL', [])
%                       CZ.n            =   noise covariance matrix of the sample mean Z.n over the periods 
%                                           size (ny+nu) x (ny+nu) x F
%                                               CZ.n(:, :) = Cov(NZ) /P
%                                               CZ.n(1:ny, 1:ny) = Cov(NY) /P 
%                                               CZ.n(ny+1:end, ny+1:end) = Cov(NU)/P 
%                                               CZ.n(1:ny, ny+1:end) = Cov(NY, NU)/P 
%                                               CZ.n(ny+1:end, 1:ny) = Cov(NU, NY)/P 
%                                           with NX the noise on X w.r.t. one period 
%
%                       CZ.NL           =   total (noise + NL dist.) covariance matrix of the sample mean Z.n over the periods 
%                                           size (ny+nu) x (ny+nu) x F
%                                               CZ.NL(:, :) = Cov(NZ) /P + Cov(ZS) 
%                                               CZ.NL(1:ny, 1:ny) = Cov(NY) /P + Cov(YS) 
%                                               CZ.NL(ny+1:end, ny+1:end) = Cov(NU) /P + Cov(US) 
%                                               CZ.NL(1:ny, ny+1:end) = Cov(NY, NU)/P + Cov(YS, US)  
%                                               CZ.NL(ny+1:end, 1:ny) = Cov(NU, NY)/P + Cov(US, YS) 
%                                           with NX the noise on X w.r.t one period, and with XS the stochastic NL distortions. 
%                                           Note: if the reference signal is not available, then  
%                                               CZ.NL(1:ny, 1:ny) = Cov(NY-G*NU) /P + Cov(YS-G*US) 
%                                               CZ.NL(ny+1:end, ny+1:end) = 0 
%
%                       CZ.S            =   covariance matrix of the stochastic nonlinear distortions in the sample mean Z.n 
%                                           over the periods  
%                                           size (ny+nu) x (ny+nu) x F
%                                               CZ.S(:, :) = Cov(ZS) 
%                                               CZ.S(1:ny, 1:ny) = Cov(YS) 
%                                               CZ.S(ny+1:end, ny+1:end) = Cov(US) 
%                                               CZ.S(1:ny, ny+1:end) = Cov(YS, US)  
%                                               CZ.S(ny+1:end, 1:ny) = Cov(US, YS)  
%                                           with XS the stochastic NL distortions 
%                                           Note: if the reference signal is not available, then  
%                                               CZ.S(1:ny, 1:ny) = Cov(YS-G*US) 
%                                               CZ.S(ny+1:end, ny+1:end) = 0 
%
%                       CZ.m_NL         =   total (noise + NL distortion) covariance matrix of the sample Z.m_NL over 
%                                           the periods and the stochastic NL distortions.   
%                                           size (ny+nu) x (ny+nu) x F
%                                               CZ.m_NL(:, :) = alpha*(Cov(NZ)/P + Cov(ZS)) 
%                                               CZ.m_NL(1:ny, 1:ny) = alpha*(Cov(NY) /P + Cov(YS)) 
%                                               CZ.m_NL(ny+1:end, ny+1:end) = alpha*(Cov(NU)/P + Cov(US)) 
%                                               CZ.m_NL(1:ny, ny+1:end) = alpha*(Cov(NY, NU)/P + Cov(YS, US))  
%                                               CZ.m_NL(ny+1:end, 1:ny) = alpha*(Cov(NU, NY)/P + Cov(US, YS)) 
%                                           with alpha < 1 the covariance reduction due to the local polynomial approximation of 
%                                           the FRM over the excited frequencies; NX the noise on X w.r.t. one period, and with 
%                                           XS the stochastic NL distortions. 
%                                           Note: if the reference signal is not available, then  
%                                               CZ.NL(1:ny, 1:ny) = alpha*(Cov(NY-G*NU) /P + Cov(YS-G*US)) 
%                                               CZ.NL(ny+1:end, ny+1:end) = 0 
%
%       Z       =   sample mean Fourier coefficients with noise and/or system transient removed, where Z stacks the output 
%                   and input spectra stacked on top of each other: Z = [Y; U]; 
%                   struct('n', [], 'm_NL', [])
%                       Z.n             =   sample mean output-input DFT spectra over the periods 
%                                           size (ny+nu) x F 
%                       Z.m_NL          =   sample mean output-input DFT spectra over the periods and the stochastic nonlinear 
%                                           distortions via the local polynomial approximation of the FRM over the excited freq.  
%
%       freq    =   frequency of the excited harmonics; size 1 x F 
%
%	    G		=	estimated frequency response matrix; size ny x nu x F
%
%       CvecG   =   covariance matrices of vec(G); size (ny*nu) x (ny*nu) x F 
%                   struct('n', [], 'NL', [], 'S')
%                       CvecG.n         =   noise covariance matrix of vec(G) 
%
%                       CvecG.NL        =   total (noise + nonlinear distortion) covariance matrix of vec(G) 
%
%                       CvecG.S         =   covariance matrix of the stochastic nonlinear distortions in vec(G) w.r.t. one
%                                           realisation and one MIMO experiment 
%
%       dof     =   actual degrees of freedom of the input-output and FRF covariance estimates = number of equivalent 
%                   independent experiments - 1 
%                   struct('n', [], 'NL', [], 'Gn', [], 'GNL', [])
%                       dof.n           =   actual degrees of freedom of the noise covariances (to reduce the uncertainty 
%                                           of the leakage suppression, the minimum dof per realisation is ny + 5)  
%                       dof.NL          =   actual degrees of freedom of the total covariance 
%                       dof.Gn          =   actual degrees of freedom of the noise covariance of vec(G) 
%                       dof.GNL         =   actual degrees of freedom of the total covariance of vec(G) 
%
%       CL      =   ± CL is the correlation length (over the frequency) of the sample mean and sample covariance 
%                   struct('n', [], 'NL', [])
%                       CL.n            =   correlation length Z.n, CZ.n, and CvecG.n over the EXCITED frequencies: 
%                                           value at frequency k is correlated with that at frequency m in [k-CL.n, k+CL.n]
%                       CL.NL           =   correlation length Z.m_NL, G, CZ.m_NL, and CvecG.NL over EXCITED frequencies 
%                                           value at frequency k is correlated with that at frequency m in [k-CL.NL, k+CL.NL]
%                   Note: CL+1 = the frequency width (in samples) of the local polynomial estimate 
%
%
%   Input
%
%       data    =   structure containing the input/output time domain signals 
%                   struct('u', [], 'y', [], 'r', [], 'N', [], 'Ts', [], 'ExcitedHarm', []) 
%                       data.u              =   input signal; size nu x N*P 
%                                               (nu inputs, P periods, N points per period) 
%                       data.y              =   ouput signal; size ny x N*P 
%                                               (ny outputs, P periods, N points per period) 
%                       data.r              =   reference signal, e.g. signal stored in the arbitrary waveform generator (optional) 
%                                               size nu x N
%                                               (nu inputs, 1 period, N points per period) 
%                                               Note: the phase spectrum varies randomly over the realisations, the experiments,
%                                                     and the frequencies, while the amplitude spectrum remains the same. 
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
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 4 November 2009 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 19 September 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = data.N;                             % number of samples per period 
[ny, NP] = size(data.y);                % ny outputs, and NP data points
nu = size(data.u, 1);                   % nu inputs
P = NP/N;                               % number of periods
fs = 1/data.Ts;                         % sampling frequency
ExcitedHarm = data.ExcitedHarm(:).';    % excited harmonic numbers
freq = ExcitedHarm*fs/N;                % frequency excited harmonics
F = length(ExcitedHarm);                % number of excited frequencies 
dof_min = ny;                           % minimal number of degrees of freedom

try
    if isempty(method)
        method = struct('order', 2, 'dof', dof_min, 'transient', 1);
    end
catch
    method = struct('order', 2, 'dof', dof_min, 'transient', 1);
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

CZ = struct('n', [], 'NL', [], 'S', [], 'm_NL', []);
CZ.NL = zeros(ny+nu, ny+nu, F);                                 % total covariance sample mean Z.n
CZ.m_NL = zeros(ny+nu, ny+nu, F);                               % sample total covariance sample mean Z.m_NL
Z = struct('n', [], 'm_NL', []);                                % output-input DFT spectra stacked on top of each other
Z.n = zeros(ny+ny, F);                                          % sample mean over the periods 
Z.m_NL = zeros(ny+ny, F);                                       % sample mean over periods and stoch. NL distortions 
G = zeros(ny, nu, F);                                           % frequency response matrix (FRM)
CvecG = struct('n', [], 'NL', [], 'S', []);
CvecG.NL = zeros(ny*nu, ny*nu, F);                              % total covariance FRM estimate
dof = struct('n', [], 'NL', [], 'Gn', [], 'GNL', []);           % degrees of freedom of the residuals used for the cov. estimates
CL = struct('n', [], 'NL', []);

if P > 1
    CZ.n = zeros(ny+nu, ny+nu, F);                              % sample noise covariance sample mean Z.n
    CZ.S = zeros(ny+nu, ny+nu, F);                              % stochastic NL contributions to Z.n
    CvecG.n = zeros(ny*nu, ny*nu, F);                           % noise covariance FRM estimate
    CvecG.S = zeros(ny*nu, ny*nu, F);                           % stochastic NL contributions to FRM estimate
end % if P > 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample means and sample noise (co-)variances output-input DFT spectra over the P periods %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample means over the periods
SelectExcited = P*ExcitedHarm + 1;                              % selection of the excited harmonics
Zfft = fft([data.y; data.u], [], 2)/(N*P);                      % scaling DTF s.t. the correct Fourier coefficients are recovered
Z.n = Zfft(:, SelectExcited);                                   % input-output Fourier coefficients

% calculate the noise (co-)variances if more than one period is available
if P > 1
        
    % calculate the sample mean and sample noise (co-)variances of the Fourier coefficients via the local polynomial approach 
    data_IO = struct('Y', [], 'ExcitedHarm', ExcitedHarm);
    method_IO = struct('order', [], 'dof', [], 'period', P, 'transient', []);
    method_IO.order = method.order;
    method_IO.dof = method.dof;               
    method_IO.transient = method.transient;
 
    % select all DFT lines (non-excited  + excited) around the excited frequency band
    StartDFT = P*ExcitedHarm(1) - (P-1);                        % at least P-1 non-excited DFT lines exist before the first harmonic
    StopDFT = P*ExcitedHarm(end) + (P-1);                       % at least P-1 non-excited DFT lines exist after the last harmonic
    SelectAll = [StartDFT:1:StopDFT] + 1; 
    data_IO.Y = Zfft(:, SelectAll);                             % noisy output-input DFT spectra      
    
    [dummy, CZm, Zm, dofn, CLn] = PeriodLocalPolyAnal(data_IO, method_IO); 
    
    Z.n = Zm;                                                   % sample mean output-input DFT spectra over the periods
    CZ.n = CZm;                                                 % sample noise covariance sample mean Z.n
    dof.n = dofn;                                               % actual degrees of freedom of the noise (co-)variance estimate
    CL.n = CLn;                                                 % correlation length Z.n and CZ.n over the frequency
    
end % P > 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRM estimate and its covariances %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data structure for the local polynomial estimate
data_FRM = struct('Y', [], 'U', [], 'freq', []);

if Reference
    % select the excited harmonics of the reference signal
    R = fft(data.r, [], 2)/N;
    R = R(:, ExcitedHarm + 1);
    data_FRM.U = R;
else % no reference signal available
    data_FRM.U = Z.n(ny+1:end,:);                               % the input DFT spectrum   
end % if reference signal available

data_FRM.Y = Z.n;
data_FRM.freq = ExcitedHarm;
if P > 1
    data_FRM.CY = CZ.n;
end % if more than one period

% method structure for the local polynomial estimate
method_FRM = struct('order', [], 'dof', []);
method_FRM.order = method.order;
method_FRM.dof = method.dof;          
if P > 1
    % the noise (and system) transient has been removed and the stochastic NL
    % distortions are periodic (=> no leakage)
    method_FRM.transient = 0;   
else
    % the noise (and system) transient has not been removed
    method_FRM.transient = 1; 
end % if more than one period
% if method.transient = 0 => no transient estimation; even if P = 1
method_FRM.transient = method_FRM.transient*method.transient;

% local polynomial estimate FRM from reference to Z and its total covariance
[CZNL, ZNL, xx, Gz, CvecGz, dofNL, CLNL] = LocalPolyAnal(data_FRM, method_FRM);

CL.NL = CLNL;                                               % correlation length Z.m_NL and CZ.m_NL over the frequency 

% FRM estimate and its total covariance
[G, CvecG] = FRF_EIV(Gz, CvecGz); 
dof.GNL = dofNL;
if P > 1
    dof.Gn = dof.n;
    % estimate stochastic nonlinear distortions on BLA
    CvecG.S = EstimStochNLdist(CvecG.NL, CvecG.n, 1);
end % if more than one period

% sample mean over periods and stoch. NL distort. and its sample total covariance
Z.m_NL = ZNL.m;
CZ.m_NL = CZNL.m;

% total covariance, and stoch. NL distort. sample mean Z.n 
if Reference
    CZ.NL = CZNL.n;                                         % total covariance output-input DFT spectra 
    dof.NL = dofNL;
    if P > 1
        CZ.S = EstimStochNLdist(CZ.NL, CZ.n, 1);
    end % if more than one period
else
    if P > 1
        % calculation output stochastic NL distortions
        CNL = CZNL.n(1:ny,1:ny,:);                          % var(NY-G*NU) + var(YS-G*US)
        Cn = zeros(ny,ny,F);                                % var(NY-G*NU)
        CS = zeros(ny,ny,F);                                % var(YS-G*US)
        for kk = 1:F
            Cn(:,:,kk) = squeeze(CZ.n(1:ny,1:ny,kk)) + squeeze(G(:,:,kk)) * squeeze(CZ.n(ny+1:end,ny+1:end,kk)) * squeeze(G(:,:,kk))' ...
                         - 2*herm(squeeze(CZ.n(1:ny,ny+1:end,kk)) * squeeze(G(:,:,kk))');
        end % kk, frequencies
        CS = EstimStochNLdist(CNL, Cn, 1);
        CZ.NL = CZNL.n;
        CZ.S(1:ny,1:ny,:) = CS;
    else
        CZ.NL(1:ny,1:ny,:) = CZNL.n(1:ny,1:ny,:);
    end % if more than one period
    dof.NL = dofNL;
end % if reference signal is available
    