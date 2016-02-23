clear all; close all; clc;
% SWEPT-SINE EX - Example of swept-sine design
%
% Experiment Parameters:
fs = 1000;              % sampling frequency   [Hz]
df = 1;                 % frequency resolution [Hz]
fl = 10;                % min excitation freq. [Hz]
fh = 400;               % max excitation freq. [Hz]
amp = 10;               % excitation amplitude [-]
[x,time,X,freq] = ssin(fs,fl,fh,df,amp);

nrofs = fs/df;
figure
subplot(211); plot(time,x(1:nrofs));
title('time domain signal'); xlabel('time [s]'); ylabel('amplitude [-]');
subplot(212); semilogx(freq,dbm(X));
title('freq domain signal'); xlabel('freq [Hz]'); ylabel('amplitude [dB]');