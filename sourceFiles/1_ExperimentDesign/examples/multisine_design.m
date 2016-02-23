%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTISINE DESIGN:
% -----------------
% Descr.:   example of multisine excitation design
%           to demonstrate multisine design factors
% Author:   Thomas Beauduin, KULeuven, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% Experiment Parameters:
fs = 1000;              % sampling frequency   [Hz]
df = 1;                 % frequency resolution [Hz]
fl = 10;                % min excitation freq. [Hz]
fh = 400;               % max excitation freq. [Hz]
amp = 10;               % excitation amplitude [-]

%% EXAMPLE 1
% full uncompressed schroeder multisine:
itp = 's';              % init phase type:  s=schroeder/r=random
ctp = 'n';              % compression type: c=comp/n=no_comp
stp = 'f';              % signal type:      f=full/ O=odd-odd
                        %                   o=odd / O2=special odd-odd
Bn=1; An=1;             % amplitude spectrum    [-]
[x1,Xs1,freqs1,Xt1,freqt1] = msin(fs,df,fl,fh,amp,itp,ctp,stp,Bn,An);

%% EXAMPLE 2
% odd compressed random mutlsine:
itp = 'r';              % init phase type:  s=schroeder/r=random
ctp = 'c';              % compression type: c=comp/n=no_comp
stp = 'o';              % signal type:      f=full/ O=odd-odd
                        %                   o=odd / O2=special odd-odd
Bn=1; An=1;             % amplitude spectrum    [-]
[x2,Xs2,freqs2,Xt2,freqt2] = msin(fs,df,fl,fh,amp,itp,ctp,stp,Bn,An);

figure
nrofs = fs/df; time=(0:1/fs:1/df-1/fs);
subplot(211); plot(time,[x1(1:nrofs),x2(1:nrofs)]);
title('time domain signal'); xlabel('time [s]'); ylabel('amplitude [-]');
subplot(212); semilogx(freqt1,[dbm(Xt1),dbm(Xt2)]);
title('freq domain signal'); xlabel('freq [Hz]'); ylabel('amplitude [dB]');

% NOTES:
% CTP: Compression type algorithm optimizes phase for 'min crest factor' 
%      to improve S/N of measurement.
% ITP: Random initial phase creates different signals in time domain
%      with identical frequency domain.
% STP: Odd signal type used for non-linear distortion analysis.
%