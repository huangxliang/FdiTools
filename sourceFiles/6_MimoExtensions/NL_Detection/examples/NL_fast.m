%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison robust and fast method for measuring the BLA, its variance,    %%
% the noise variance, and the variance of the stochastic nonlinear          %%
% distortions from noisy input/output measurements using random phase       %%
% multisines                                                                %%
%                                                                           %%
% Illustrated on a discrete-time Wiener-Hammerstein system consisting       %%
% of the cascade of                                                         %%
%                                                                           %%
%   1. a second order bandpass system:  G1(z^-1)                            %%
%   2. a static nonlinear function:     alpha*x^2 + beta*x^3                %%
%   3. a second order bandpass system:  G2(z^-1)                            %%
%                                                                           %%
% Note: this example needs the "DesignMultisineExcitations" toolbox         %%
%                                                                           %%
% Rik Pintelon                                                              %% 
% December 3, 2007                                                          %%
%                                                                           %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition random phase multisine excitation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 700e6;                       % frequency first excited harmonic
fmax = 1100e6;                      % frequency last excited harmonic
fres = 1e6;                         % frequency resolution excited harmonics
fs = 4e9;                           % sampling frequency
frat = 1.02;                        % ratio consecutive excited harmonics; only for logarithmic spacing
Nblock = 3;                         % number of consecutive harmonics out of which one is randomly eliminated

% linear spacing excited harmonics
Spacing = 'lin';  

% logarithmic spacing of the excited harmonics
% Spacing = 'log'; fmin = 10e6;

% odd excited harmonics only
TypeMultisine = 'odd';

% odd and even excited harmonics
% TypeMultisine = 'full';

if strcmp(TypeMultisine, 'full')
    fres = fres/2;
end % if


DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fs = fs;
DefFreq.fres = fres;
if strcmp(Spacing, 'log')
    DefFreq.frat = frat;                % only needed for logarithmic spacing
end % if logarithmic spacing

[ExcitedHarm, N, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, TypeMultisine);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation measurement %
%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 20;                                     % number of realisations random phase multisine
% M = 1;
P = 6;                                      % number of consecutive periods

fstart = 600e6;                             % first measured frequency
fstop = 1200e6;                             % last measured frequency
% fstart = 1e6;
% fstop = 1200e6*2;                         % last measured frequency

freqall = ([0:1:N-1]*fs/N).';               % all DFT frequencies over full unit circle
qall = exp(-sqrt(-1)*2*pi*freqall/fs);      % z^-1 over full unit circle
qOdE = exp(-sqrt(-1)*2*pi*ExcitedHarm/N);   % z^-1 at odd excited frequencies
if Spacing == 'log'
    if strcmp(TypeMultisine, 'odd')
        fstart = fres / 2;
    else % full multisine
        fstart = fres ;
    end % if
    fstop = fmax * 3;
end % if
if strcmp(TypeMultisine, 'odd')
    freqmeas = ([fstart:fres/2:fstop]).';       % all measured frequencies
else % full multisine
    freqmeas = ([fstart:fres:fstop]).';       % all measured frequencies
end % if
MeasHarm = freqmeas/(fs/N);                 % measured harmonics
Fall = length(MeasHarm);                    % number of measured frequencies
Uall = zeros(M, P, Fall);                   % all input spectra at the measured frequencies
Yall = zeros(M, P, Fall);                   % all output spectra at the measured frequencies


% discrete-time LTI1
a1 = [1, -0.2, 0.9];
b1 = [1];
G1all = polyval(fliplr(b1), qall)./polyval(fliplr(a1), qall);
G1 = polyval(fliplr(b1), qOdE)./polyval(fliplr(a1), qOdE);

% polynomial static nonlinearity
alfa = 0.1;                                % coefficient  x^2
beta = 0.001;                              % coefficient x^3
if Spacing == 'log'
    beta = 0.1;
end % if

% discrete-time LTI2
a2 = [1, -0.5, 0.9];
b2 = [0, 1, 0.5];
G2all = polyval(fliplr(b2), qall)./polyval(fliplr(a2), qall);
G2 = polyval(fliplr(b2), qOdE)./polyval(fliplr(a2), qOdE);

% input/output noise standard deviations
sigmau = 0.1;      % input noise std
sigmay = 0.2;      % output noise std

for ii=1:M
    
    % generation one period random phase multisine with rms value one 
    % and with the same harmonic content
    u = CalcMultisine(ExcitedHarm, N);
    
    % response to the first LTI system
    U = fft(u)/sqrt(N);
    X = U.*G1all;
    x = real(ifft(X)*sqrt(N));
    
    % response of the static nonlinearity
    z = x + alfa*x.^2 + beta*x.^3;
    
    % response to the second LTI system
    Z = fft(z)/sqrt(N);
    Y = Z.*G2all;
    
    % collect measurements
    Uall(ii,:,:) = repmat((U(MeasHarm+1)).', P, 1);
    Yall(ii,:,:) = repmat((Y(MeasHarm+1)).', P, 1);
   
end

% noiseless input
U0all = Uall;
% add white input and output noise to the measurements
Uall = Uall + (randn(size(Uall)) + sqrt(-1)*randn(size(Uall)))* sigmau/sqrt(2);
Yall = Yall + (randn(size(Yall)) + sqrt(-1)*randn(size(Yall))) * sigmay/sqrt(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required variables:                                                    %
%   ExcitedHarm     =   excited odd harmonics in multisine               %
%                           size: F x 1                                  %
%   MeasHarm        =   all measured harmonics                           %
%                           size: Fall x 1                               %
%   N               =   number of samples in one period of the multisine %
%   fs              =   sampling frequency                               %
%   Uall, Yall      =   M x P X Fall input output spectra where          %
%                           M = number of realisations multisine         %
%                           P = number of consecutive periods multisine  %
%                           Fall = number of all measured frequencies    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of all measured harmonics
MeasHarm = round(freqmeas/(fs/N));

% NL analysis on IO spectra
[Y, Yc, U, G, freq] = Fast_NL_Anal(Yall, Uall, ExcitedHarm, MeasHarm, fs, N);

% plot results
RealNum = 5;            % number of realisation to be plotted
FigNum = 1;             % number of first figure
Plot_Fast_NL_Anal(Y, Yc, U, G, freq, Spacing, RealNum, FigNum);
