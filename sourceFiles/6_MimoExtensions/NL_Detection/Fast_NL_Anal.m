function [Y, Yc, U, G, freq] = Fast_NL_Anal(Yall, Uall, ExcitedHarm, MeasHarm, fs, N);
%
%       function [Y, Yc, U, G, freq] = Fast_NL_Anal(Yall, Uall, ExcitedHarm, MeasHarm, fs, N);
%       
%       Usage:
%
%           Quantification and qualification (odd and/or even) of the nonlinear distortions in FRF measurements using 
%           full or odd random phase multisines with random harmonic grid. These signals can be generated using the 
%           "DesignMultisineExcitations" toolbox. The analysis starts from P >= 2 periods of the steady state response to 
%           a full or odd random phase multisine with random harmonic grid. M > 1 experiments with independent random phase 
%           realisations of the multisines can be handled if the multisines have the same rms value, the same amplitude spectrum, 
%           and the same harmonic content. 
%
%       References:
%
%           Pintelon R., G. Vandersteen, L. De Locht, Y. Rolain and J. Schoukens (2004). Experimental characterization of 
%           operational amplifiers: a system identification approach - Part I: theory and simulations, IEEE Trans. Instrum. Meas., 
%           vol. 53, no. 3, pp. 854-862.
%
%           Schoukens, J., J. Lataire, R. Pintelon, G. Vandersteen, and T. Dobrowiecki (2009). Robustness issues of the equivalent 
%           linear representation of a nonlinear system, IEEE Trans. Instrum. and Meas., vol. 58, no. 5, pp. 1737-1745.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA). 
%
%
%   OUTPUT
%
%       Y       =   struct{'mean', 'std'} containing mean value and standard deviation output spectrum
%                   Y.mean  =   struct{'E', 'NE'} containing the mean value of the output spectrum over P consecutive periods; size M x P x NumberOfFrequencies
%                               Y.mean.E    =   mean output spectrum at the excited (odd) harmonic frequencies
%                               Y.mean.NE   =   struct{'even', 'odd'} for ODD multisines
%                                               struct{'all', 'inband', outband'} for FULL multisines 
%                                               Y.mean.NE.even      =   struct {'all', 'inband', 'outband'}
%                                                                       Y.mean.NE.even.all      =   mean output spectrum at all the non-excited even harmonics
%                                                                       Y.mean.NE.even.inband   =   mean output spectrum at the inband non-excited even harmonics
%                                                                       Y.mean.NE.even.outband  =   mean output spectrum at the outband non-excited even harmonics
%                                               Y.mean.NE.odd       =   struct {'all', 'inband', 'outband'}
%                                                                       Y.mean.NE.odd.all       =   mean output spectrum at all the non-excited odd harmonics
%                                                                       Y.mean.NE.odd.inband    =   mean output spectrum at the inband non-excited odd harmonics
%                                                                       Y.mean.NE.odd.outband   =   mean output spectrum at the outband non-excited odd harmonics
%                                               Y.mean.NE.all       =   mean output spectrum at all non-excited harmonics 
%                                               Y.mean.NE.inband    =   mean output spectrum at the inband non-excited harmonics 
%                                               Y.mean.NE.outband	=   mean output spectrum at the outband non-excited harmonics 
%                   Y.std   =   struct{'E', 'NE'} containing the standard deviation of the mean value of the output spectrum over P consecutive periods
%                               Y.std.E     =   std output spectrum at the excited (odd) harmonic frequencies
%                               Y.std.NE    =   struct{'even', 'odd'} for ODD multisines 
%                                               struct{'all', 'inband', outband'} for FULL multisines 
%                                               Y.std.NE.even       =   struct {'all', 'inband', 'outband'}
%                                                                       Y.std.NE.even.all       =   std output spectrum at all the non-excited even harmonics
%                                                                       Y.std.NE.even.inband    =   std output spectrum at the inband non-excited even harmonics
%                                                                       Y.std.NE.even.outband   =   std output spectrum at the outband non-excited even harmonics
%                                               Y.std.NE.odd        =   struct {'all', 'inband', 'outband'}
%                                                                       Y.std.NE.odd.all        =   std output spectrum at all the non-excited odd harmonics
%                                                                       Y.std.NE.odd.inband     =   std output spectrum at the inband non-excited odd harmonics
%                                                                       Y.std.NE.odd.outband    =   std output spectrum at the outband non-excited odd harmonics
%                                               Y.std.NE.all        =   std output spectrum at all non-excited harmonics 
%                                               Y.std.NE.inband     =   std output spectrum at the inband non-excited harmonics 
%                                               Y.std.NE.outband	=   std output spectrum at the outband non-excited harmonics 
%
%       Yc      =   structure {'mean', 'std'} containing the inband mean value and standard deviation corrected output spectrum; size M x P x NumberOfFrequencies
%                   corrected for the presence of energy at non-excited input harmonics
%                   Yc.mean =   struct{'E', 'NE'} containing the inband mean corrected output spectrum
%                               Yc.mean.E   =   mean corrected output spectrum at the excited (odd) harmonics (same as uncorrected)
%                               Yc.mean.NE  =   struct{'even', 'odd'} for ODD multisines 
%                                               Yc.mean.NE.even =   mean corrected output spectrum at non-excited inband even harmonics
%                                               Yc.mean.NE.odd  =   mean corrected output spectrum at non-excited inband odd harmonics
%                                               mean corrected output spectrum at the non-excited harmonics for FULL multisines 
%                   Yc.std  =   struct{'E', 'NE'} containing the inband standard deviation of the corrected output spectrum
%                               Yc.std.E    =   std corrected output spectrum at excited odd harmonics (same as uncorrected)
%                               Yc.std.NE   =   struct{'even', 'odd'} for ODD multisines 
%                                               Yc.std.NE.even  =   std corrected output spectrum at non-excited inband even harmonics
%                                               Yc.std.NE.odd   =   std corrected output spectrum at non-excited inband odd harmonics
%                                               std corrected output spectrum at the non-excited harmonics for FULL multisines 
%
%       U       =   structure {'mean', 'std'} containing mean value and standard deviation input spectrum; size M x P x F
%                   U.mean  =   struct{'E', 'NE'} containing the mean value of the input spectrum over P consecutive periods
%                               U.mean.E    =   mean input spectrum at excited odd harmonic frequencies
%                               U.mean.NE   =   struct{'even', 'odd'} for ODD multisines 
%                                               struct{'all', 'inband', outband'} for FULL multisines 
%                                               U.mean.NE.even      =   struct {'all', 'inband', 'outband'}
%                                                                       U.mean.NE.even.all      =   mean input spectrum at all the non-excited even harmonics
%                                                                       U.mean.NE.even.inband   =   mean input spectrum at the inband non-excited even harmonics
%                                                                       U.mean.NE.even.outband  =   mean input spectrum at the outband non-excited even harmonics
%                                               U.mean.NE.odd       =   struct {'all', 'inband', 'outband'}
%                                                                       U.mean.NE.odd.all       =   mean input spectrum at all the non-excited odd harmonics
%                                                                       U.mean.NE.odd.inband    =   mean input spectrum at the inband non-excited odd harmonics
%                                                                       U.mean.NE.odd.outband   =   mean input spectrum at the outband non-excited odd harmonics
%                                               U.mean.NE.all       =   mean input spectrum at all non-excited harmonics 
%                                               U.mean.NE.inband    =   mean input spectrum at the inband non-excited harmonics 
%                                               U.mean.NE.outband	=   mean input spectrum at the outband non-excited harmonics 
%                   U.std   =   struct{'E', 'NE'} containing the standard deviation of the mean value of the input spectrum over P consecutive periods
%                               U.std.E     =   std input spectrum at excited odd harmonic frequencies
%                               U.std.NE    =   struct{'even', 'odd'} for ODD multisines 
%                                               struct{'all', 'inband', outband'} for FULL multisines 
%                                               U.std.NE.even       =   struct {'all', 'inband', 'outband'}
%                                                                       U.std.NE.even.all       =   std input spectrum at all the non-excited even harmonics
%                                                                       U.std.NE.even.inband    =   std input spectrum at the inband non-excited even harmonics
%                                                                       U.std.NE.even.outband   =   std input spectrum at the outband non-excited even harmonics
%                                               U.std.NE.odd        =   struct {'all', 'inband', 'outband'}
%                                                                       U.std.NE.odd.all        =   std input spectrum at all the non-excited odd harmonics
%                                                                       U.std.NE.odd.inband     =   std input spectrum at the inband non-excited odd harmonics
%                                                                       U.std.NE.odd.outband    =   std input spectrum at the outband non-excited odd harmonics
%                                               U.std.NE.all        =   std input spectrum at all non-excited harmonics 
%                                               U.std.NE.inband     =   std input spectrum at the inband non-excited harmonics 
%                                               U.std.NE.outband	=   std input spectrum at the outband non-excited harmonics 
%
%       G       =   structure {'all', 'mean', 'stdn', 'stdNL', 'stds'} containing mean value and standard deviations FRF at excited odd harmonics
%                   G.all   =   Frequency Response Function (FRF) for all realisations and periods; size M x P x F
%                   G.mean  =   mean value FRF over the P consecutive periods
%                   G.stdn  =   struct{'E', 'NE'} containing the noise standard deviation of the mean FRF value over the P consecutive periods
%                               G.stdn.E    =   noise std calculated from the excited frequencies
%                               G.stdn.NE   =   noise std calculated from the non-excited odd frequencies (= extrapolation)
%                               Note: a difference between G.stdn.E and G.stdn.NE indicates a non-stationary behaviour
%                   G.stdNL =   total standard deviation mean value FRF over P consecutive periods:
%                               G.stdNL.^2  =   noise variance + stochatic NL distortions
%                   G.stds =    standard deviation stochastic NL distortions
%                               G.stds  =   sqrt(|G.stdNL.^2 - G.stdn.^2|)
%                               Note: G.stds is a measure of the order of magnitude of the bias contribution of the NL distortions to the
%                                     best linear approximation
%
%       freq    =   struct {'E', 'NE'} containing the excited and non-excited frequencies
%                   freq.E  =   vector containing the excited (odd) harmonic frequencies
%                   freq.NE =   struct{'even', 'odd'} for ODD multisines
%                               struct{'all', 'inband', outband'} for FULL multisines 
%                               freq.NE.even        =   struct{'all', 'inband', outband'}
%                                                       freq.NE.even.all        =   all non-excited even harmonics
%                                                       freq.NE.even.inband     =   inband non-excited even harmonics
%                                                       freq.NE.even.outband    =   outband non-excited even harmonics
%                               freq.NE.odd         =   struct{'all', 'inband', outband'}
%                                                       freq.NE.odd.all         =   all non-excited odd harmonics
%                                                       freq.NE.odd.inband      =   inband non-excited odd harmonics
%                                                       freq.NE.odd.outband     =   outband non-excited odd harmonics
%                               freq.NE.all         =   all non-excited harmonics 
%                               freq.NE.inband      =   inband non-excited harmonics 
%                               freq.NE.outband     =   outband non-excited harmonics 
%
%   INPUT
%
%       Yall        =   M x P x Fmeas (or P x Fmeas) matrix containing the output spectrum
%                       M = number of independent random phase realisations with same rms value, amplitude spectrum, and harmonic grid
%                       P = number of consecutive periods
%                       Fmeas = number of measured frequencies
%
%       Uall        =   M x P x Fmeas (or P x Fmeas) matrix containing the input spectrum
%                       M = number of independent random phase realisations with same rms value, amplitude spectrum, and harmonic grid
%                       P = number of consecutive periods
%                       Fmeas = number of measured frequencies
%
%       ExcitedHarm =   excited (odd) harmonics; size Fexcited x 1
%       MeasHarm    =   measured harmonics; size Fmeas x 1     
%       fs          =   sampling frequency
%       N           =   number of points in one period
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, March 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 5 December 2007
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation variables and structures OUTPUT parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExcitedHarm = ExcitedHarm(:);
% test whether only the odd harmonics are excited
OddMultisine = isempty(find((ExcitedHarm - floor(ExcitedHarm/2)*2) == 0, 1));

MeasHarm = MeasHarm(:);
TheSize = size(Yall);
if length(TheSize) == 3
    M = TheSize(1);
    P = TheSize(2);
    F = TheSize(3);
else
    M = 1;
    P = TheSize(1);
    F = TheSize(2);
    dummy = zeros(1, P, F);
    dummy(1, :, :) = Yall;
    Yall = dummy;
    dummy(1, :, :) = Uall;
    Uall = dummy;
end % if

Y = struct('mean', [], 'std', []);
Y.mean = struct('E', [], 'NE', []);
if OddMultisine
    Y.mean.NE = struct('even', [], 'odd', []);
    Y.mean.NE.even = struct('all', [], 'inband', [], 'outband', []);
    Y.mean.NE.odd = Y.mean.NE.even;
else % full multisine
    Y.mean.NE = struct('all', [], 'inband', [], 'outband', []);
end % if
Y.std = Y.mean;
freq = Y.mean;
U = Y;
Yc = struct('mean', [], 'std', []);
Yc.mean = struct('E', [], 'NE', []);
if OddMultisine
    Yc.mean.NE = struct('even', [], 'odd', []);
end % if odd multisine
Yc.std = Yc.mean;
G = struct('all', [], 'mean', [], 'stdn', [], 'stdNL', [], 'stds', []);
G.stdn = struct('E', [], 'NE', []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% harmonic content and frequencies %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% harmonic content
NonExcitedHarm = HarmonicContent(MeasHarm, ExcitedHarm);

% frequencies
freq.E = ExcitedHarm * fs/N;
if OddMultisine
    freq.NE.even.all = NonExcitedHarm.even.all * fs/N;
    freq.NE.even.inband = NonExcitedHarm.even.inband * fs/N;
    freq.NE.even.outband = NonExcitedHarm.even.outband * fs/N;
    freq.NE.odd.all = NonExcitedHarm.odd.all * fs/N;
    freq.NE.odd.inband = NonExcitedHarm.odd.inband * fs/N;
    freq.NE.odd.outband = NonExcitedHarm.odd.outband * fs/N;
else % full multisine
    freq.NE.all = NonExcitedHarm.all * fs/N;    
    freq.NE.inband = NonExcitedHarm.inband * fs/N;    
    freq.NE.outband = NonExcitedHarm.outband * fs/N;    
end % if


%%%%%%%%%%%%%%%%%%%%%%%%
% input-output spectra %
%%%%%%%%%%%%%%%%%%%%%%%%

Uallm = mean(Uall, 2);                                          % mean over consecutive periods
Yallm = mean(Yall, 2);                                          % mean over consecutive periods
resUall = Uall - repmat(Uallm, [1, P, 1]);
resYall = Yall - repmat(Yallm, [1, P, 1]);
stdUallm = squeeze((sum(abs(resUall.^2), 2)/(P-1)/P).^0.5);     % std mean value
stdYallm = squeeze((sum(abs(resYall.^2), 2)/(P-1)/P).^0.5);     % std mean value
varYU = squeeze(sum(resYall.*conj(resUall), 2)/(P-1)/P);        % covariance mean values
Uallm = squeeze(Uallm);
Yallm = squeeze(Yallm);
if M == 1
    Uallm = Uallm.';
    stdUallm = stdUallm.';
    Yallm = Yallm.';
    stdYallm = stdYallm.';
    varYU = varYU.';
end % if

% covariance at the excited frequencies;
% to be used for calculating the variance of G over consecutive periods
varYU = varYU(:, ExcitedHarm - MeasHarm(1) + 1);

% input spectrum, mean and standard deviation
U.mean.E = Uallm(:, ExcitedHarm - MeasHarm(1) + 1);
if OddMultisine
    U.mean.NE.even.all = Uallm(:, NonExcitedHarm.even.all - MeasHarm(1) + 1);
    U.mean.NE.even.inband = Uallm(:, NonExcitedHarm.even.inband - MeasHarm(1) + 1);
    U.mean.NE.even.outband = Uallm(:, NonExcitedHarm.even.outband - MeasHarm(1) + 1);
    U.mean.NE.odd.all = Uallm(:, NonExcitedHarm.odd.all - MeasHarm(1) + 1);
    U.mean.NE.odd.inband = Uallm(:, NonExcitedHarm.odd.inband - MeasHarm(1) + 1);
    U.mean.NE.odd.outband = Uallm(:, NonExcitedHarm.odd.outband - MeasHarm(1) + 1);
else % full multisine
    U.mean.NE.all = Uallm(:, NonExcitedHarm.all - MeasHarm(1) + 1);
    U.mean.NE.inband = Uallm(:, NonExcitedHarm.inband - MeasHarm(1) + 1);
    U.mean.NE.outband = Uallm(:, NonExcitedHarm.outband - MeasHarm(1) + 1);
end % if OddMultisine

U.std.E = stdUallm(:, ExcitedHarm - MeasHarm(1) + 1);
if OddMultisine
    U.std.NE.even.all = stdUallm(:, NonExcitedHarm.even.all - MeasHarm(1) + 1);
    U.std.NE.even.inband = stdUallm(:, NonExcitedHarm.even.inband - MeasHarm(1) + 1);
    U.std.NE.even.outband = stdUallm(:, NonExcitedHarm.even.outband - MeasHarm(1) + 1);
    U.std.NE.odd.all = stdUallm(:, NonExcitedHarm.odd.all - MeasHarm(1) + 1);
    U.std.NE.odd.inband = stdUallm(:, NonExcitedHarm.odd.inband - MeasHarm(1) + 1);
    U.std.NE.odd.outband = stdUallm(:, NonExcitedHarm.odd.outband - MeasHarm(1) + 1);
else % full multisine
    U.std.NE.all = stdUallm(:, NonExcitedHarm.all - MeasHarm(1) + 1);
    U.std.NE.inband = stdUallm(:, NonExcitedHarm.inband - MeasHarm(1) + 1);
    U.std.NE.outband = stdUallm(:, NonExcitedHarm.outband - MeasHarm(1) + 1);
end % if OddMultisine

% output spectrum, mean and standard deviation
Y.mean.E = Yallm(:, ExcitedHarm - MeasHarm(1) + 1);
if OddMultisine
    Y.mean.NE.even.all = Yallm(:, NonExcitedHarm.even.all - MeasHarm(1) + 1);
    Y.mean.NE.even.inband = Yallm(:, NonExcitedHarm.even.inband - MeasHarm(1) + 1);
    Y.mean.NE.even.outband = Yallm(:, NonExcitedHarm.even.outband - MeasHarm(1) + 1);
    Y.mean.NE.odd.all = Yallm(:, NonExcitedHarm.odd.all - MeasHarm(1) + 1);
    Y.mean.NE.odd.inband = Yallm(:, NonExcitedHarm.odd.inband - MeasHarm(1) + 1);
    Y.mean.NE.odd.outband = Yallm(:, NonExcitedHarm.odd.outband - MeasHarm(1) + 1);
else % full multisine
    Y.mean.NE.all = Yallm(:, NonExcitedHarm.all - MeasHarm(1) + 1);
    Y.mean.NE.inband = Yallm(:, NonExcitedHarm.inband - MeasHarm(1) + 1);
    Y.mean.NE.outband = Yallm(:, NonExcitedHarm.outband - MeasHarm(1) + 1);
end % if OddMultisine

Y.std.E = stdYallm(:, ExcitedHarm - MeasHarm(1) + 1);
if OddMultisine
    Y.std.NE.even.all = stdYallm(:, NonExcitedHarm.even.all - MeasHarm(1) + 1);
    Y.std.NE.even.inband = stdYallm(:, NonExcitedHarm.even.inband - MeasHarm(1) + 1);
    Y.std.NE.even.outband = stdYallm(:, NonExcitedHarm.even.outband - MeasHarm(1) + 1);
    Y.std.NE.odd.all = stdYallm(:, NonExcitedHarm.odd.all - MeasHarm(1) + 1);
    Y.std.NE.odd.inband = stdYallm(:, NonExcitedHarm.odd.inband - MeasHarm(1) + 1);
    Y.std.NE.odd.outband = stdYallm(:, NonExcitedHarm.odd.outband - MeasHarm(1) + 1);
else % full multisine
    Y.std.NE.all = stdYallm(:, NonExcitedHarm.all - MeasHarm(1) + 1);
    Y.std.NE.inband = stdYallm(:, NonExcitedHarm.inband - MeasHarm(1) + 1);
    Y.std.NE.outband = stdYallm(:, NonExcitedHarm.outband - MeasHarm(1) + 1);
end % if OddMultisine


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Response Function (FRF)                                          %
% 1. mean value over the P consecutive periods                               %
% 2. noise standard deviation mean value over the P consecutive periods      %
% 3. prediction total distortion on FRF from level non-excited odd harmonics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G.all = Yall(:, :, ExcitedHarm - MeasHarm(1) + 1) ./ Uall(:, :, ExcitedHarm - MeasHarm(1) + 1);

% mean value and noise standard deviation FRF over P consecutive periods at the excited odd harmonics
% the mean FRF is calculated as the ratio of the mean values of the spectra in order to reduce bias
% introduced by input noise in case of poor input SNR
Gallm = Y.mean.E ./ U.mean.E;
stdGallm = abs(Gallm) .* ((Y.std.E.^2)./abs(Y.mean.E.^2) + (U.std.E.^2)./abs(U.mean.E.^2) ...
                           - 2*real(varYU./(Y.mean.E .* conj(U.mean.E)))).^0.5;

G.mean = Gallm;
G.stdn.E = stdGallm;            % noise std from excited odd harmonics

% calculation FRF at inband non-excited harmonics via linear interpolation
if OddMultisine
    NonExcitedInBandHarm = sort([NonExcitedHarm.even.inband; NonExcitedHarm.odd.inband]);
else % full multisine
    NonExcitedInBandHarm = NonExcitedHarm.inband;
end % if OddMultisine
GallNEib = interp1(ExcitedHarm, Gallm.', NonExcitedInBandHarm).';       % linear interpolation is performed over the columns
                                                                        % GallOdNE contains NaN where interpolation does not exist

% correction output spectrum at inband non-excited harmonics for presence distortions at input spectrum
Ycall = Yall;
TheIndex = NonExcitedInBandHarm - MeasHarm(1) + 1;
Ycall(:, :, TheIndex) = Ycall(:, :, TheIndex) - permute(repmat(GallNEib, [1 1 P]) .* permute(Uall(:, :, TheIndex), [1 3 2]), [1 3 2]);

% mean and standard deviation corrected output spectrum
Ycallm = squeeze(mean(Ycall, 2));                   % mean value over periods
stdYcallm = squeeze(std(Ycall, 0, 2))/sqrt(P);      % std mean value over periods
if M == 1
    Ycallm = Ycallm.';
    stdYcallm = stdYcallm.';
end % if

Yc.mean.E = Ycallm(:, ExcitedHarm - MeasHarm(1) + 1);
if OddMultisine
    Yc.mean.NE.even = Ycallm(:, NonExcitedHarm.even.inband - MeasHarm(1) + 1);
    Yc.mean.NE.odd = Ycallm(:, NonExcitedHarm.odd.inband - MeasHarm(1) + 1);
else % full multisine
    Yc.mean.NE = Ycallm(:, NonExcitedHarm.inband - MeasHarm(1) + 1);
end % if OddMultisine

Yc.std.E = stdYcallm(:, ExcitedHarm - MeasHarm(1) + 1);
if OddMultisine
    Yc.std.NE.even = stdYcallm(:, NonExcitedHarm.even.inband - MeasHarm(1) + 1);
    Yc.std.NE.odd = stdYcallm(:, NonExcitedHarm.odd.inband - MeasHarm(1) + 1);
else % full multisine
    Yc.std.NE = stdYcallm(:, NonExcitedHarm.inband - MeasHarm(1) + 1);
end % if OddMultisine

% estimate total distortion on FRF measurement via extrapolation level distortion at inband non-excited (odd) harmonics to excited odd harmonics
% note: take absolute value !!
if OddMultisine
    Ysall = interp1(NonExcitedHarm.odd.inband, abs(Yc.mean.NE.odd).', ExcitedHarm).';       % linear interpolation is performed for each column
else % full multisine
    Ysall = interp1(NonExcitedHarm.inband, abs(Yc.mean.NE).', ExcitedHarm).';               % linear interpolation is performed for each column
end % if OddMultisine
G.stdNL = Ysall ./ abs(U.mean.E);

% estimate standard deviation stochastic NL distortion
G.stds = (abs(G.stdNL.^2 - G.stdn.E.^2)).^0.5;

% estimate noise level on FRF measurement via extrapolation noise level at inband non-excited harmonics to excited (odd) harmonics
if OddMultisine
    AllInBandNE = [NonExcitedHarm.odd.inband; NonExcitedHarm.even.inband];
    AllYcstdNE = [Yc.std.NE.odd.'; Yc.std.NE.even.'];       % row index is the frequency, column index the number of realisations
    [AllInBandNE, TheIndex] = sort(AllInBandNE);
    AllYcstdNE = AllYcstdNE(TheIndex, :);
else % full multisine
    AllInBandNE = NonExcitedHarm.inband;
    AllYcstdNE = Yc.std.NE.';                               % row index is the frequency, column index the number of realisations
end % if OddMultisine
Ynall = interp1(AllInBandNE, AllYcstdNE, ExcitedHarm).';    % linear interpolation is performed for each column
G.stdn.NE = Ynall ./ abs(U.mean.E);
