clear all; close all; clc;
% NONPARA_EX - Nonparmetric Model Identification example.
% Example Data:
load('exdata_ol.mat');

% UNFINISHED

% STEP 1
% data pretreatment and visual inspection:
x=pretreat(iq_ref,nrofs,fs);
[y,time]=pretreat(pos_t,nrofs,fs);

% STEP 2
% estimator selection and frf model identification:
[X,Y,FRFs,FRFn,freq] = time2frf_ml(x,y,fs,fl,fh,nrofs);
[FRFh,freq2,coh1,X2,Y2] = time2frf_h1(x,y,fs,fl,fh,fr,1);
[FRFl,freq3,coh2] = time2frf_log(x,y,fs,fl,fh,fr,hanning(nrofs),nrofs*0.1);

figure
subplot(211)
semilogx(freq,[dbm(FRFs),dbm(FRFn)]);
title('FRF estimate'); ylabel('magnitude [dB]');
subplot(212)
semilogx(freq,phs(FRFs))
xlabel('freq [Hz]'); ylabel('phase [deg]');

figure
subplot(211)
semilogx(freq3,[dbm(FRFh),dbm(FRFl)]);
title('FRF estimate'); ylabel('magnitude [dB]');
subplot(212)
semilogx(freq3,[phsw(FRFh),phsw(FRFl)])
xlabel('freq [Hz]'); ylabel('phase [deg]');