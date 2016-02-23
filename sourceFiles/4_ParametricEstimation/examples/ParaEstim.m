clear all; close all; clc;
% PARAESTIM - parametric model estimation example.
% Nonparametric Estimation - ML:
load('exdata_ol.mat');
x=pretreat(iq_ref,nrofs,fs); y=pretreat(pos_t,nrofs,fs);
[X,Y,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,y,fs,fl,fh,nrofs);

% Parametric Estimation:
n=4;                        % model order of denominator polynomial
mh=3; ml=0;                 % model orders of numerator polynomial
r=(1:440);                  % frequency lines with information
relvar=1e-10;               % relative variation of costfunction (stop)
iter=1e3;                   % maximum number of iterations (stop)
GN = 0;                     % optimization type (Levenberg-Marquardt)
cORd = 'c';                 % identifaction model type (continuous)
FRF_W = ones(size(FRFs));   % least squares weigting function
relax = 1;                  % relaxation factor for btls

[Bn_nls,An_nls,Bn_ls,An_ls,Bn_wls,An_wls] = ...
nllsfdi(FRFs(r),freq(r),FRF_W(r),n,mh,ml,iter,relvar,GN,cORd,fs);
SYSnls = tf(Bn_nls,An_nls); FRFnls=squeeze(freqresp(SYSnls,freq*2*pi));
SYSls = tf(Bn_ls,An_ls); FRFls=squeeze(freqresp(SYSls,freq*2*pi));
SYSwls = tf(Bn_wls,An_wls); FRFwls=squeeze(freqresp(SYSwls,freq*2*pi));

[Bn_ml,An_ml] =...
mlfdi(X(r),Y(r),freq(r),sX2(r),sY2(r),cXY(r),n,mh,ml,iter,relvar,cORd,fs);
SYSml = tf(Bn_ml,An_ml); FRFml=squeeze(freqresp(SYSml,freq*2*pi));

[Bn_btls,An_btls,Bn_gtls,An_gtls] = ...
btlsfdi(Y(r),X(r),freq(r),n,mh,ml,sY2(r),sX2(r),cXY(r),iter,relax,relvar);
SYSbtls = tf(Bn_btls,An_btls); FRFbtls=squeeze(freqresp(SYSbtls,freq*2*pi));
SYSgtls = tf(Bn_gtls,An_gtls); FRFgtls=squeeze(freqresp(SYSgtls,freq*2*pi));

figure
subplot(211)
semilogx(freq,[dbm(FRFs),dbm(FRFls),dbm(FRFnls),dbm(FRFml)...
              ,dbm(FRFbtls),dbm(FRFgtls)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on
legend('FRF','LSE','NLS','MLE','BTLS','GTLS')
subplot(212)
semilogx(freq,[phs(FRFs),phs(FRFls),phs(FRFnls),phs(FRFml)...
              ,phs(FRFbtls),phs(FRFgtls)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
