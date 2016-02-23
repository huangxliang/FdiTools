clear all; close all; clc;
% BINARY EX - Example of PBRS design
% Experiment Parameters:

% UNFINNISHED

fs = 1000;              % sampling frequency   [Hz]
df = 1;                 % frequency resolution [Hz]
fl = 10;                % min excitation freq. [Hz]
fh = 400;               % max excitation freq. [Hz]
amp = 10;               % excitation amplitude [-]
log2N = 30;
bitno = fs/df;
startnum = [];
[x,nextstnum]=prbs(log2N,bitno,startnum);

figure
nrofs = fh/df; time=(0:1/fs:1/df-1/fs);
X=t2f(x,nrofs); freq=fs*(0:1:nrofs/2-1)'/nrofs;
subplot(211); plot(x(1:nrofs));
title('time domain signal'); xlabel('time [s]'); ylabel('amplitude [-]');
subplot(212); semilogx(freq,dbm(X));
title('freq domain signal'); xlabel('freq [Hz]'); ylabel('amplitude [dB]');

x=-ones(nrofs,1);
x(mod([1:nrofs].^2,nrofs)+1)=1;
x(1)= 1;
figure
X=t2f(x,nrofs); freq=fs*(0:1:nrofs/2-1)'/nrofs;
subplot(211); plot(x(1:nrofs));
title('time domain signal'); xlabel('time [s]'); ylabel('amplitude [-]');
subplot(212); semilogx(freq,dbm(X));
title('freq domain signal'); xlabel('freq [Hz]'); ylabel('amplitude [dB]');
