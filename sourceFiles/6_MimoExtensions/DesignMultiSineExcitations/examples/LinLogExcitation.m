%
% Calculation of random phase multisine excitations with random harmonic grid:
%
%   1. uniform distribution harmonics
%   2. logarithmic distribution harmonics
%
% for the following types of multisines
%
%   1. odd multisines (the even harmonics are not excited)
%   2. full multisines (all harmonics are excited)
%
% Rik Pintelon, 4 October 2011
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition odd lin tone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 700e6;               % start frequency multisine in Hz
fmax = 1100e6;              % stop frequency multisine in Hz
fres = 1e6;                 % desired frequency resolution in Hz
fs = 4e9;                   % sampling frequency in Hz
Nblock = 3;                 % one out of three consecutive harmonics is randomly eliminated
Spacing = 'lin';            % linear frequency spacing
MultiType = 'odd';          % no even excited harmonics
% fmin = fres/2;          % lowpass design

DefFreq.fs = fs;
DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fres = fres;

% calculation of the random harmonic grid of the odd multisine with linear spacing excited harmonics 
[ExcitedHarmOddLin, N_odd_lin, NewDefFreqOddLin] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);
length(ExcitedHarmOddLin)

% one realisation odd random phase multisine with a given random harmonic grid
% this routine is used each time a random phase realisation is needed (the harmonic 
% grid remains the same) 
OddLinSignal = CalcMultisine(ExcitedHarmOddLin, N_odd_lin);

% plot DFT spectrum
OddLinSpec = fft(OddLinSignal)/sqrt(N_odd_lin);
SelectAll = [1:N_odd_lin/2+1].';
freqAll = (SelectAll-1)*fs/N_odd_lin;
OddLinSpec = OddLinSpec(SelectAll);
figure(1)
subplot(211)
plot(OddLinSignal);
xlabel('samples')
ylabel('time signal')
title('Odd lin tone')
subplot(212)
plot(freqAll, db(OddLinSpec),'+')
xlabel('Frequency (Hz)')
ylabel('DFT spectrum')
zoom on
shg


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition odd log tone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 1e2;                 % start frequency multisine in Hz
fmax = 1e4;                 % stop frequency multisine in Hz
fres = 10;                  % desired frequency resolution in Hz
frat = 1.015;               % ratio frequencies of two consecutive excited harmonics
fs = 1.2e5;                 % sampling frequency in Hz
Nblock = 3;                 % one out of three consecutive harmonics is randomly eliminated
Spacing = 'LOG';            % logarithmic spacing excited harmonics
MultiType = 'odd';          % no even excited harmonics
% fmin = fres/2;          % lowpass design

DefFreq.fs = fs;
DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fres = fres;
DefFreq.frat = frat;

% calculation of the random harmonic grid of the odd multisine with logarithmic spacing excited harmonics 
[ExcitedHarmOddLog, N_odd_log, NewDefFreqOddLog] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);
length(ExcitedHarmOddLog)

% one realisation odd random phase multisine with a given random harmonic grid
% this routine is used each time a random phase realisation is needed (the harmonic 
% grid remains the same) 
OddLogSignal = CalcMultisine(ExcitedHarmOddLog, N_odd_log);

% plot DFT spectrum
OddLogSpec = fft(OddLogSignal)/sqrt(N_odd_log);
SelectAll = [1:N_odd_log/2+1].';
freqAll = (SelectAll-1)*fs/N_odd_log;
OddLogSpec = OddLogSpec(SelectAll);
figure(2)
subplot(211)
plot(OddLogSignal);
xlabel('samples')
ylabel('time signal')
title('Odd log tone')
subplot(212)
semilogx(freqAll, db(OddLogSpec),'+')
xlabel('Frequency (Hz)')
ylabel('DFT spectrum')
zoom on
shg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition full lin tone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 700e6;               % start frequency multisine in Hz
fmax = 1100e6;              % stop frequency multisine in Hz
fres = 1e6;                 % desired frequency resolution in Hz
fs = 4e9;                   % sampling frequency in Hz
Nblock = 3;                 % one out of three consecutive harmonics is randomly eliminated
Spacing = 'lin';            % linear spacing excited harmonics
MultiType = 'full';         % odd and even harmonics are excited
% fmin = fres/2;          % lowpass design

DefFreq.fs = fs;
DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fres = fres;

% calculation of the random harmonic grid of the full multisine with linear distribution excited harmonics 
[ExcitedHarmFullLin, N_full_lin, NewDefFreqFullLin] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);
length(ExcitedHarmFullLin)

% one realisation odd random phase multisine with a given random harmonic grid
% this routine is used each time a random phase realisation is needed (the harmonic 
% grid remains the same) 
FullLinSignal = CalcMultisine(ExcitedHarmFullLin, N_full_lin);

% plot DFT spectrum
FullLinSpec = fft(FullLinSignal)/sqrt(N_full_lin);
SelectAll = [1:N_full_lin/2+1].';
freqAll = (SelectAll-1)*fs/N_full_lin;
FullLinSpec = FullLinSpec(SelectAll);
figure(3)
subplot(211)
plot(FullLinSignal);
xlabel('samples')
ylabel('time signal')
title('Full lin tone')
subplot(212)
plot(freqAll, db(FullLinSpec),'+')
xlabel('Frequency (Hz)')
ylabel('DFT spectrum')
zoom on
shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition full log tone %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmin = 1e2;                 % start frequency multisine in Hz
fmax = 1e4;                 % stop frequency multisine in Hz
fres = 10;                  % desired frequency resolution in Hz
frat = 1.015;               % ratio frequencies of two consecutive excited harmonics
fs = 1.2e5;                 % sampling frequency in Hz
Nblock = 3;                 % one out of three consecutive harmonics is randomly eliminated
Spacing = 'LOG';            % logarithmic spacing excited harmonics
MultiType = 'full';         % odd and even harmonics are excited
% fmin = fres/2;          % lowpass design

DefFreq.fs = fs;
DefFreq.fmin = fmin;
DefFreq.fmax = fmax;
DefFreq.fres = fres;
DefFreq.frat = frat;

% calculation of the random harmonic grid of the full multisine with logarithmic distribution excited harmonics 
[ExcitedHarmFullLog, N_full_log, NewDefFreqFullLog] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);
length(ExcitedHarmFullLog)

% one realisation odd random phase multisine with a given random harmonic grid 
% this routine is used each time a random phase realisation is needed (the harmonic 
% grid remains the same) 
FullLogSignal = CalcMultisine(ExcitedHarmFullLog, N_full_log);

% plot DFT spectrum
FullLogSpec = fft(FullLogSignal)/sqrt(N_full_log);
SelectAll = [1:N_full_log/2+1].';
freqAll = (SelectAll-1)*fs/N_full_log;
FullLogSpec = FullLogSpec(SelectAll);
figure(4)
subplot(211)
plot(FullLogSignal);
xlabel('samples')
ylabel('time signal')
title('Full log tone')
subplot(212)
semilogx(freqAll, db(FullLogSpec),'+')
xlabel('Frequency (Hz)')
ylabel('DFT spectrum')
zoom on
shg
