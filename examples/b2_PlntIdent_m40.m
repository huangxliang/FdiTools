%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4th ORDER PLANT IDENTIFICATION:
% -------------------------------
% Descr.:   Axial 1e mode transfer functions (Yuna)
% System:   general ball-screw-driven machine-tool axis
% Author:   Thomas Beauduin, University of Tokyo, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
load('data\hres\DhresL_m40pR.mat');

% STEP 1: DATA
% data pretreatment and visual inspection
x=pretreat(iq_ad,nrofs,fs,1);
ym=pretreat(theta_m,nrofs,fs,1);
[yl,time]=pretreat((pos_t.*1e-6),nrofs,fs,1);
figure, r0=(1:nrofs); m=4; rm=(nrofs*m+1:(m+1)*nrofs);
plotyy(time(r0),x(rm),time(r0),[ym(rm),yl(rm)]); 
legend('Torque','Angle','Position')


% STEP 2: FRF
% frf estimator selection and identification:
[Xm,Ym,FRFsm,FRFnm,freq,sX2m,sY2m,cXYm,sCRm] = time2frf_ml(x,ym,fs,fl,fh,nrofs);
[Xl,Yl,FRFsl,FRFnl,freq,sX2l,sY2l,cXYl,sCRl] = time2frf_ml(x,yl,fs,fl,fh,nrofs);
figure
subplot(221), semilogx(freq,[dbm(FRFsm),dbm(FRFnm),dbm(sCRm)]);
title('Motor-side'); ylabel('magnitude [dB]'); xlim([fl,fh]);
subplot(223), semilogx(freq,phs(FRFsm))
xlabel('freq [Hz]'); ylabel('phase [deg]'); xlim([fl,fh]);
subplot(222), semilogx(freq,[dbm(FRFsl),dbm(FRFnl),dbm(sCRl)]);
title('Load-side'); ylabel('magnitude [dB]'); xlim([fl,fh]);
subplot(224), semilogx(freq,phs(FRFsl))
xlabel('freq [Hz]'); ylabel('phase [deg]'); xlim([fl,fh]);

% STEP 3 - MOTOR
% parametric estimators comparison
[X,Y,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,ym,fs,fl,fh,nrofs);
n=5;                        % model order of denominator polynomial
mh=3; ml=0;                 % model orders of numerator polynomial
r=(1:280);                  % frequency lines with information
relvar=1e-10;               % relative variation of costfunction (stop)
iter=1e3;                   % maximum number of iterations (stop)
GN = 0;                     % Levenberg-Marquardt optimization
cORd = 'c';                 % continuous model identifaction
FRF_W = ones(size(FRFs));   % least squares weigting function
relax = 1;                  % relaxation factor for btls
[Bn_nls,An_nls,Bn_ls,An_ls] = ...
nlsfdi(FRFs(r),freq(r),FRF_W(r),n,mh,ml,iter,relvar,GN,cORd,fs);
SYSnls = tf(Bn_nls,An_nls); FRFnls=squeeze(freqresp(SYSnls,freq*2*pi));
SYSls = tf(Bn_ls,An_ls); FRFls=squeeze(freqresp(SYSls,freq*2*pi));
[Bn_ml,An_ml] =...
mlfdi(X(r),Y(r),freq(r),sX2(r),sY2(r),cXY(r),n,mh,ml,iter,relvar,cORd,fs);
SYSml = tf(Bn_ml,An_ml); FRFml=squeeze(freqresp(SYSml,freq*2*pi));
[Bn_btls,An_btls,Bn_gtls,An_gtls] = ...
btlsfdi(Y(r),X(r),freq(r),n,mh,ml,sY2(r),sX2(r),cXY(r),iter,relax,relvar);
SYSbtls = tf(Bn_btls,An_btls); FRFbtls=squeeze(freqresp(SYSbtls,freq*2*pi));
SYSgtls = tf(Bn_gtls,An_gtls); FRFgtls=squeeze(freqresp(SYSgtls,freq*2*pi));
figure, subplot(221)
semilogx(freq,[dbm(FRFs),dbm(FRFls),dbm(FRFnls)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on, legend('FRF','WLS','NLS')
title('Deterministic Methods')
subplot(223)
semilogx(freq,[phs(FRFs),phs(FRFls),phs(FRFnls)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
subplot(222)
semilogx(freq,dbm(FRFs),freq,dbm(FRFml),'c',freq,dbm(FRFbtls),'m','LineWidth',2)
ylabel('Amplitude [dB]'), grid on, legend('FRF','MLE','BTLS')
title('Stochastic Methods')
subplot(224)
semilogx(freq,phs(FRFs),freq,phs(FRFml),'c',freq,phs(FRFbtls),'m','LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
SYS4 = SYSml; FRF4 = FRFml;

% STEP 4 - MOTOR
% fine tune system identification based on requirements
[B2,A2]=mlfdi(X,Y,freq,sX2,sY2,cXY,2,0,0,iter,relvar,cORd,fs);
SYS2=tf(B2,A2); FRF2=squeeze(freqresp(SYS2,freq*2*pi));
[~,pp,kp]=zpkdata(SYS2,'v');
pp=sort(pp); pp(2)=-0; pp(1)=real(pp(1));
SYSm2=zpk([],pp,kp); FRFm2=squeeze(freqresp(SYSm2,freq*2*pi));

[zp,pp,~]=zpkdata(SYS4,'v');
for k=1:mh;
    if real(zp(k))>0; zp(k)=-real(zp(k))+1i*imag(zp(k)); end
end
zp=zp(2:end);
pp=sort(pp); pp(1)=-0; pp(2)=-0; pp=pp(1:4);
SYSm4 = zpk(zp,pp,kp); FRFm4=squeeze(freqresp(SYSm4,freq*2*pi));
figure;
subplot(221); semilogx(freq,[dbm(FRFs),dbm(FRF2),dbm(FRFm2),dbm(FRFn)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on, title('model order 2')
subplot(223); semilogx(freq,[phs(FRFs),phs(FRFm2)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
subplot(222); semilogx(freq,[dbm(FRFs),dbm(FRF4),dbm(FRFm4),dbm(FRFn)],'LineWidth',2)
grid on, legend('signal','estim','model','noise'), title('model order 4')
subplot(224); semilogx(freq,[phs(FRFs),phs(FRF4),phs(FRFm4)],'LineWidth',2)
xlabel('frequency [Hz]'),grid on

% STEP 3 - LOAD
% parametric estimators comparison
[X,Y,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,yl,fs,fl,fh,nrofs);
n=6;                        % model order of denominator polynomial
mh=3; ml=0;                 % model orders of numerator polynomial
r=(1:290);                  % frequency lines with information
relvar=1e-10;               % relative variation of costfunction (stop)
iter=5e2;                   % maximum number of iterations (stop)
GN = 0;                     % Levenberg-Marquardt optimization
cORd = 'c';                 % continuous model identifaction
FRF_W = ones(size(FRFs));   % least squares weigting function
relax = 1;                  % relaxation factor for btls
[Bn_nls,An_nls,Bn_ls,An_ls] = ...
nlsfdi(FRFs(r),freq(r),FRF_W(r),n,mh,ml,iter,relvar,GN,cORd,fs);
SYSnls = tf(Bn_nls,An_nls); FRFnls=squeeze(freqresp(SYSnls,freq*2*pi));
SYSls = tf(Bn_ls,An_ls); FRFls=squeeze(freqresp(SYSls,freq*2*pi));
[Bn_ml,An_ml] =...
mlfdi(X(r),Y(r),freq(r),sX2(r),sY2(r),cXY(r),n,mh,ml,iter,relvar,cORd,fs);
SYSml = tf(Bn_ml,An_ml); FRFml=squeeze(freqresp(SYSml,freq*2*pi));
[Bn_btls,An_btls,Bn_gtls,An_gtls] = ...
btlsfdi(Y(r),X(r),freq(r),n,mh,ml,sY2(r),sX2(r),cXY(r),iter,relax,relvar);
SYSbtls = tf(Bn_btls,An_btls); FRFbtls=squeeze(freqresp(SYSbtls,freq*2*pi));
SYSgtls = tf(Bn_gtls,An_gtls); FRFgtls=squeeze(freqresp(SYSgtls,freq*2*pi));
figure, subplot(221)
semilogx(freq,[dbm(FRFs),dbm(FRFls),dbm(FRFnls)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on, legend('FRF','WLS','NLS')
title('Deterministic Estimation Methods')
subplot(223)
semilogx(freq,[phs(FRFs),phs(FRFls),phs(FRFnls)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
subplot(222)
semilogx(freq,dbm(FRFs),freq,dbm(FRFml),'c',freq,dbm(FRFbtls),'m','LineWidth',2)
ylabel('Amplitude [dB]'), grid on, legend('FRF','MLE','BTLS')
title('Stochastic Estimation Methods')
subplot(224)
semilogx(freq,phs(FRFs),freq,phs(FRFml),'c',freq,phs(FRFbtls),'m','LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
SYS6 = SYSml; FRF6 = FRFml;

% STEP 4 - LOAD
% fine tune system identification based on requirements
r = (1:50);
[B2,A2]=mlfdi(X(r),Y(r),freq(r),sX2(r),sY2(r),cXY(r),2,0,0,iter,relvar,cORd,fs);
SYS2=tf(B2,A2); FRF2=squeeze(freqresp(SYS2,freq*2*pi));
[~,pp,kp]=zpkdata(SYS2,'v');
pp=sort(pp); pp(2)=-0; pp(1)=real(pp(1));
SYSl2=zpk([],pp,kp); FRFl2=squeeze(freqresp(SYSl2,freq*2*pi));

[zp,pp,~]=zpkdata(SYS6,'v');
for k=1:mh;
    if real(zp(k))>0; zp(k)=-real(zp(k))+1i*imag(zp(k)); end
end
zp=zp(2:end); kp=0.4175e+08*(12/2/pi);
pp=sort(pp); pp(1)=-0; pp(2)=-0; 
SYSl4 = zpk([],pp(1:4),kp); FRFl4=squeeze(freqresp(SYSl4,freq*2*pi));
SYSl6 = zpk(zp,pp(1:6),kp*1.25); FRFl6=squeeze(freqresp(SYSl6,freq*2*pi));
figure;
subplot(231); semilogx(freq,[dbm(FRFs),dbm(FRF2),dbm(FRFl2),dbm(FRFn)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on, title('model order 2')
subplot(234); semilogx(freq,[phs(FRFs),phs(FRFl2)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
subplot(232); semilogx(freq,[dbm(FRFs),dbm(FRF6),dbm(FRFl4),dbm(FRFn)],'LineWidth',2)
grid on, title('model order 4')
subplot(235); semilogx(freq,[phs(FRFs),phs(FRF6),phs(FRFl4)],'LineWidth',2)
xlabel('frequency [Hz]'),grid on
subplot(233); semilogx(freq,[dbm(FRFs),dbm(FRF6),dbm(FRFl6),dbm(FRFn)],'LineWidth',2)
grid on, title('model order 6'), legend('signal','estim','model','noise')
subplot(236); semilogx(freq,[phs(FRFs),phs(FRF6),phs(FRFl6)],'LineWidth',2)
xlabel('frequency [Hz]'),grid on

% STEP 5 - COMBINE: N2
% combine the design to SIMO and display
[zm,pm,km]=zpkdata(SYSm2,'v'); sys1 = SYSm2;
[zl,pl,kl]=zpkdata(SYSl2,'v'); sys2 = zpk(zl,pm,kl);
sysn2 = ss([sys1;sys2]);
SSm2=squeeze(freqresp(sysn2(1),freq*2*pi));
SSl2=squeeze(freqresp(sysn2(2),freq*2*pi));
figure;
subplot(221); semilogx(freq,[dbm(FRFsm),dbm(FRFm2),dbm(SSm2)],'LineWidth',2)
subplot(223); semilogx(freq,[phsw(FRFsm),phsw(FRFm2),phsw(SSm2)],'LineWidth',2)
subplot(222); semilogx(freq,[dbm(FRFsl),dbm(FRFl2),dbm(SSl2)],'LineWidth',2)
subplot(224); semilogx(freq,[phs(FRFsl),phs(FRFl2),phs(SSl2)],'LineWidth',2)

% STEP 5 - COMBINE: N4
% combine the design to SIMO and display
[zl,pl,kl]=zpkdata(SYSl4,'v');
[zm,pm,km]=zpkdata(SYSm4,'v');
da=real(zm(1,1))/real(pm(3,1)); db=imag(zm(2,1))/imag(pm(3,1));
zm(1,1)=( real(pl(3,1))*da )+1i*( imag(pl(3,1))*db );
zm(2,1)=( real(pl(3,1))*da )-1i*( imag(pl(3,1))*db );
sys1 = zpk(zm,pl,km); sys2 = zpk(zl,pl,kl);
sysn4 = ss([sys1;sys2]);
SSm4=squeeze(freqresp(sysn4(1),freq*2*pi));
SSl4=squeeze(freqresp(sysn4(2),freq*2*pi));
figure;
subplot(221); semilogx(freq,[dbm(FRFsm),dbm(FRFm4),dbm(SSm4)],'LineWidth',2)
subplot(223); semilogx(freq,[phsw(FRFsm),phsw(FRFm4),phsw(SSm4)],'LineWidth',2)
subplot(222); semilogx(freq,[dbm(FRFsl),dbm(FRFl4),dbm(SSl4)],'LineWidth',2)
subplot(224); semilogx(freq,[phs(FRFsl),phs(FRFl4),phs(SSl4)],'LineWidth',2)
break
% STEP 6 - SAVE
% save results & create plant model header for csim
P=load('data\PsysData.mat');
P.sys40n2 = sysn2; P.sys40n4=sysn4;
P.sys40n4t=SYSl4; P.sys40n4m=SYSm4; P.sys40n6t=SYSl6;
P.FRFs40=[FRFsm,FRFsl]; P.FRFn40=[FRFnm,FRFnl]; P.freq40 = freq;
save('data\PsysData.mat','-struct','P')

file = 'firm\A_StageFloat\data\system_sim_par.h';
title = 'BALL-SCREW SETUP MODEL';
phdr = 'Open-Loop Table Position:';
sysd = c2d(P.sys40n4(2),1/fs,'zoh');
[S.Aplnt,S.Bplnt,S.Cplnt,S.Dplnt]=ssdata(sysd);
mdl2hdr(file,title,{phdr},{S});

% TF to SIMO:
% design SIMO state-space by combining 2TF's
% choose the poles from motor-side, beter designed
% in function of the zeros.

% PRINCIPAL COORDINATES:
% use the modal data to create principal coordinates
% state space formulations; to get state-space with
% physical meaning as required

% SS2SS, CANON:
% only use the matlab tools ss2ss, canon to get standard
% forms of the state-space representation, not arbirary ones.

% NUMERICALLY STABLE:
% scale the state-space matrix with balance or prescale
% to get numerically stable ss for embedded system run

% GRAMMIAN SYS REDUCTION:
% during selection and simplification of system, consider 
% using grammian reduction in state-space to avoid gain fitting
