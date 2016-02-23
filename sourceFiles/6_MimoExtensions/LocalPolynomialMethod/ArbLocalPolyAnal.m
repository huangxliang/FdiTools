function [CZ, Z, freq, G, CvecG, dof, CL] = ArbLocalPolyAnal(data, method);
%
%		Using noisy input-output data of 1 MIMO experiment with arbitrary excitations on an nu
%		input and ny output system, the following quatities are calculated:
%
%           (i)     The generalized sample means and sample covariances (noise or noise + nonlinear distortions)
%                   of the input/output DFT spectra with optimal noise and/or system transient suppression. 
%
%           (ii)	The nu x ny frequency response matrix (FRM) and its covariance (noise or noise + nonlinear distortions). 
%       
%       Notes:
%
%           (i)     To handle the noisy input, noisy output case (errors-in-variables framework) a known  
%                   reference signal should be available. If not (data.r is empty), then the estimates are biased due to the input 
%                   noise. If the input is known exactly (generalized output error framework), then no reference signal is needed. 
%
%           (ii)    A local polynomial approximation of the FRM and the noise and system transient term is made. 
%
%           (iii)   The input signal may be exactly zero in (parts of) the frequency band of interest.
%
%           (iv)	The algorithm can handle data sets from different experiments obtained under the same operational conditions.
%
%           (v)     For nonlinear systems the FRM is the best linear approximation and the covariance is the 
%                   sum of the noise covariance and the covariance of the stochastic nonlinear distortions. 
%
%           (vi)    If no input data is provided, then the algorithm simplifies to nonparametric time series 
%                   analysis (noise power spectrum estimation).
%
%       Warning:    The estimated frequency response matrix is meaningless in those frequency bands were 
%                   the input is exactly zero.
%
%       References:
%
%                   Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and 
%                   FRF models for multivariable systems - Part I: theory, Mechanical Systems and Signal Processing,
%                   vol. 24, no. 3, pp. 573-595.
%
%                   Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and 
%                   FRF models for multivariable systems - Part II: extensions, applications, Mechanical Systems and 
%                   Signal Processing, vol. 24, no. 3, pp. 596-616.
%
%                   Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%                   IEEE Press-Wiley, Piscataway (USA). 
%       
%
%	function [CZ, Z, freq, G, CvecG, dof, CL] = ArbLocalPolyAnal(data, method);
%
%
%	Output parameters
%
%		CZ      =	struct('n', [], 'm', [], 'm_nt', [])
%                       CZ.n	:   sample covariance matrix output-input noise (+ nonlinear distortions); 
%                                   size (ny+nu) x (ny+nu) x F
%                                       Usage: calculation uncertainty bounds on the FRM 
%                       CZ.m	:	sample covariance matrix sample mean output-input; size (ny+nu) x (ny+nu) x F
%                                       Usage: weighting in frequency domain maximum likelihood estimation 
%                       CZ.m_nt	:	sample covariance matrix sample mean output-input with transient removed; 
%                                   size (ny+nu) x (ny+nu) x F
%                                       Usage: weighting in frequency domain maximum likelihood estimation 
%                       Remark  :   for output error problems (the input is known exactly) only the 
%                                   ny x ny upper blocks are nonzero 
%
%		Z      =   struct('m', [], 'm_nt', [])
%                       Z.m     :   sample mean output-input Z = [Y; U]; size (ny+nu) x F
%                                       Usage: frequency domain maximum likelihood estimation 
%                       Z.m_nt  :   sample mean Z with transient removed, size (ny+nu) x F
%                                       Usage: - frequency domain maximum likelihood estimation
%                                              - leakage free output DFT spectra 
%
%       freq    =   frequency vector; size 1 x F 
%
%		TZ      =	sum plant and noise transient contribution in Z; size (ny+nu) x F
%                       Usage: - transient removal in sample mean Zm(k) 
%                              - transient removal in output spectrum Z(k)
%                       Remark  :   for output error problems (the input is known exactly) only the 
%                                   first ny rows are nonzero 
%
%       G       =   estimated frequency response matrix; size ny x nu x F
%
%       CvecG   =   covariance matrix vec(G) (noise + nonlinear distortions); size (ny*nu) x (ny*nu) x F 
%
%       dof     =   actual degrees of freedom of the residuals used for the covariance estimate 
%                   (= equivalent number of independent experiments-1) 
%
%       CL      =   ± CL is the correlation length (over the frequency) of the sample mean and sample covariance:
%                   Y(k) is correlated with Y(m) with m in [k-CL, k+CL]
%                   Note: CL+1 = the frequency width (in samples) of the local polynomial estimate 
%
%
%   Input
%
%       data    =   structure containing the input/output time domain signals 
%                   struct('u', [], 'y', [], 'N', [], 'Ts', []) 
%                       data.u              =   (concatenated) input signals; size nu x N 
%                                               (nu inputs, P periods, N data points) 
%                       data.y              =   (concatenated) ouput signals; size ny x N 
%                                               (ny outputs, P periods, N data points) 
%                       data.r              =   (concatenated) reference signals (optional), e.g. signal stored in the
%                                               arbitrary waveform generator; size nu x N 
%                                               	Usage: - if data.r is empty, then data.u should be exactly known 
%                                                            (= generalized output error problem) 
%                                                          - if data.r is given, then data.u might be noisy 
%                                                            (= errors-in-variables problem) 
%                       data.N              =   array containing the number of time domain samples in each 
%                                               data record (optional; default data.N = N)
%                                               size 1 X number of concatenated data records (experiments) 
%                       data.Ts             =   sampling period
%
%       method  =   structure containing the parameters of the method used (optional) 
%                   struct('order', [], 'dof', [], 'transient', [], 'startfreq', [], 'stopfreq', [], 'step', []) 
%                       method.order        =   order of the local polynomial approximation (optional; default 2) 
%                       method.dof          =   degrees of freedom of the (co-)variance estimates = equivalent number of 
%                                               independent experiments - 1 (optional; default ny) 
%                       method.transient    = 	determines the estimation of the transient term(s) (optional; default 1)  
%                                                   1: transient term is estimated 
%                                                   0: no transient term is estimated 
%                       method.startfreq    =   start frequency in Hz of the output parameters (optional; default fs/N); 
%                                               the frequency is rounded towards the smallest DFT frequency larger than startfreq: 
%                                               startfreq = fs/N*ceil(startfreq/(fs/N))  
%                                               
%                       method.stopfreq     =   stop frequency in Hz of the output parameters (optional; default (N/2-1)*fs/N) 
%                                               the frequency is rounded towards the largest DFT frequency smaller than stopfreq: 
%                                               stopfreq = fs/N*floor(stopfreq/(fs/N))  
%                                                
%                       method.step         =   integer number that determines the frequency resolution at which  
%                                               the output parameters are calculated (optional; default 1):
%                                                   the maximal frequency resolution is fs/N (step = 1) 
%                                                   otherwise the frequency resolution is step*fs/N 
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 15 September 2011 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 12 October 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny, N] = size(data.y);                 % ny outputs, and N data points
nu = size(data.u, 1);                   % nu inputs

% number of samples in each concatenated data record 
try
	if isempty(data.N)
		data.N = N;
	end
catch
	data.N = N;
end % try
Nconcat = length(data.N);               % number of concatenated data records

fs = 1/data.Ts;                         % sampling frequency
fres = fs/N;                            % frequency resolution FRF measurement
dof_min = ny;                           % minimal number of degrees of freedom

try
    if isempty(method)
        method = struct('order', 2, 'dof', dof_min, 'transient', 1, 'step', 1);
        % startfreq and stopfreq are normalised on the frequency resolution
    end
catch
    method = struct('order', 2, 'dof', dof_min, 'transient', 1, 'step', 1);
    % startfreq and stopfreq are normalised on the frequency resolution
end % try

% order local polynomial approximation
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

% the start frequency is normalised on the frequency resolution
% (value should be larger than DC)
try
	if isempty(method.startfreq)
		method.startfreq = 1;                   % normalised on the frequency resolution
    else
        method.startfreq = ceil(method.startfreq/fres);
	end
catch
	method.startfreq = 1;
end % try
if (method.startfreq < 1) || (method.startfreq > N/2-1)
    method.startfreq = 1;
end % if

% the stop frequency is normalised on the frequency resolution 
% value should be smaller than Nyquist 
try
	if isempty(method.stopfreq)
		method.stopfreq = N/2-1;                   % normalised on the frequency resolution
    else
        method.stopfreq = floor(method.stopfreq/fres);
	end
catch
	method.stopfreq = N/2-1;
end % try
if (method.stopfreq > N/2-1) || (method.stopfreq < 1)
    method.stopfreq = N/2-1;
end % if
if method.stopfreq < method.startfreq
    method.stopfreq = method.startfreq;
end % if

% integer number determining the frequency resolution of the calculations
try
	if isempty(method.step)
		method.step = 1; 
	end
catch
	method.step = 1;
end % try
if method.step < 0
    method.step = 1;
end % if

try
    if isempty(data.r)
        Reference = 0;                  % generalized output error; input is known, noisy output
    else
        Reference = 1;                  % errors-in-variables; noisy input and output
    end % if
catch
    data.r = [];
    Reference = 0;                      % generalized output error; input is known, noisy output
end % try

% selection of the DFT lines of the output parameters
SelectDFTlines = [method.startfreq:method.step:method.stopfreq]; 

% frequency vector
freq = SelectDFTlines*fs/N;             % frequency selected DFT lines
F = length(freq);                       % number of selected DFT lines 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional inputs in case of concatenated data records %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uadd = zeros(Nconcat-1, N);             % empty in case of 1 data record
PositionImpulse = 1;                    % position impulse in additional input 
for ii = 1:Nconcat-1
    PositionImpulse = PositionImpulse + data.N(ii);
    uadd(ii, PositionImpulse) = 1;
end % ii


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFT spectra at the requested DFT frequencies %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = fft(data.u, [], 2)/sqrt(N);
U = U(:, SelectDFTlines+1);                 % column position = DFT line number + 1
clear data.u

Y = fft(data.y, [], 2)/sqrt(N);
Y = Y(:, SelectDFTlines+1);                 % column position = DFT line number + 1
clear data.y

if Reference
    R = fft(data.r, [], 2)/sqrt(N);
    R = R(:, SelectDFTlines+1);
    clear data.r
end % if reference signal available

if Nconcat > 1
    if method.transient == 1
        Uadd = fft(uadd, [], 2)/sqrt(N);
        Uadd = Uadd(:, SelectDFTlines+1);
    else
        % if no transients are estimated then the concatenated data sets
        % are handled as one data set
        Uadd = [];
        Nconcat = 1;
    end % if
    clear uadd
else
    Uadd = [];
end % if concatenated data records


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation FRF via the LocalPolyAnal function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.freq = freq;               % frequencies at which the analysis is performed
data.nu_add = Nconcat-1;        % indicates that the data is the result of the concatenation of Nconcat data sets
method.step = 1;                % FRF is calculated at all the requested frequencies

if Reference == 1   % errors-in-variables
    
    data.U = [R; Uadd];
    data.Y = [Y; U];
    
    % generalized sample mean and sample covariance Z = [Y; U]
    % and estimate FRF from R to Z and its covariance
    [CZ, Z, TZ, GRZ, CvecGRZ, dof, CL] = LocalPolyAnal(data, method); 
    
    % FRM from U to Y and its covariance
    [G, CvecG] = FRF_EIV(GRZ, CvecGRZ);
       
else                % generalized output error
    
    data.U = [U; Uadd];
    data.Y = Y;
    
    % generalized sample mean and sample covariance Y
    % and estimate FRF from U to Y and its covariance
    [CY, Y, TY, G, CvecG, dof, CL] = LocalPolyAnal(data, method); 
    
    % generalized sample mean of Z = [Y; U]
    Z = struct('m', [], 'm_nt', []);
    Z.m = [Y.m; U];
    Z.m_nt = [Y.m_nt; U];
    
    % generalized sample covariance of Z = [Y; U]
    CZ = struct('n', zeros(ny+nu, ny+nu, F), 'm', zeros(ny+nu, ny+nu, F), 'm_nt', zeros(ny+nu, ny+nu, F));
    CZ.n(1:ny, 1:ny, :) = CY.n;
    CZ.m(1:ny, 1:ny, :) = CY.m;
    CZ.m_nt(1:ny, 1:ny, :) = CY.m_nt;
    
end % if Reference

