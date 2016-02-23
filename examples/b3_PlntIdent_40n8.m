%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4th ORDER PLANT IDENTIFICATION:
% -------------------------------
% Descr.:   Axial 1e mode transfer functions (Yuna)
% System:   general ball-screw-driven machine-tool axis
% Author:   Thomas Beauduin, University of Tokyo, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
L=load('data\hres\Dhres_40n8l.mat');
H=load('data\hres\Dhres_40n8h.mat');

% STEP 1 - ALL DATA
% data pretreatment and visual inspection:
[xl,tl]=pretreat(L.iq_ad,L.nrofs,L.fs,1); [xh,th]=pretreat(H.iq_ad,H.nrofs,H.fs,1);
yl=pretreat(L.pos_t.*1e-6,L.nrofs,L.fs,1); yh=pretreat(H.pos_t.*1e-6,H.nrofs,H.fs,1);
yl2=pretreat(L.theta_m,L.nrofs,L.fs,1); yh2=pretreat(H.theta_m,H.nrofs,H.fs,1);
figure; 
subplot(221),plot(tl,[xl,yl]); legend('input','output')
subplot(222),plot(th,[xh,yh]); legend('input','output')
subplot(223),plot(tl,[xl,yl2]); legend('input','output')
subplot(224),plot(th,[xh,yh2]); legend('input','output')

% STEP 2 - TABLE
% estimator selection and frf model identification:
[Xl,Yl,FRFsl,FRFnl,freql,sX2l,sY2l,cXYl,sCRl] = time2frf_ml(xl,yl,L.fs,L.fl,L.fh,L.nrofs);
[Xh,Yh,FRFsh,FRFnh,freqh,sX2h,sY2h,cXYh,sCRh] = time2frf_ml(xh,yh,H.fs,H.fl,H.fh,H.nrofs);
fs = H.fs; fl=L.fl; fh=H.fh; nrofs=L.nrofs+H.nrofs;
X=vertcat(Xl,Xh); Y=vertcat(Yl,Yh); FRFs=vertcat(FRFsl,FRFsh);
FRFn=vertcat(FRFnl,FRFnh); freq=vertcat(freql,freqh); sX2=vertcat(sX2l,sX2h); 
sY2=vertcat(sY2l,sY2h); cXY=vertcat(cXYl,cXYh);  sCR=vertcat(sCRl,sCRh);
figure, subplot(211), semilogx(freq,[dbm(FRFs),dbm(FRFn)]);
title('FRF estimate'); ylabel('magnitude [dB]');
subplot(212), semilogx(freq,phs(FRFs));
xlabel('freq [Hz]'); ylabel('phase [deg]');

% STEP 3 - TABLE
% parametric estimation
n=8;                        % model order of denominator polynomial
mh=3; ml=0;                 % model orders of numerator polynomial
r=(1:390);                  % frequency lines with information
relvar=1e-10;               % relative variation of costfunction (stop)
iter=4e2;                   % maximum number of iterations (stop)
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
figure
subplot(211), semilogx(freq,[dbm(FRFs),dbm(FRFls),dbm(FRFnls),...
               dbm(FRFml),dbm(FRFbtls),dbm(FRFgtls)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on, legend('FRF','LSE','NLS','MLE','BTLS','GTLS')
subplot(212), semilogx(freq,[phs(FRFs),phs(FRFls),phs(FRFnls),phs(FRFml)...
              ,phs(FRFbtls),phs(FRFgtls)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
FRF4t=FRFs; FRF4tn = FRFn; SYSt = SYSbtls;

% STEP 2 - MOTOR
% estimator selection and frf model identification:
[Xl,Yl,FRFsl,FRFnl,freql,sX2l,sY2l,cXYl,sCRl] = time2frf_ml(xl,yl2,L.fs,L.fl,L.fh,L.nrofs);
[Xh,Yh,FRFsh,FRFnh,freqh,sX2h,sY2h,cXYh,sCRh] = time2frf_ml(xh,yh2,H.fs,H.fl,H.fh,H.nrofs);
fs = H.fs; fl=L.fl; fh=H.fh; nrofs=L.nrofs+H.nrofs;
X=vertcat(Xl,Xh); Y=vertcat(Yl,Yh); FRFs=vertcat(FRFsl,FRFsh);
FRFn=vertcat(FRFnl,FRFnh); freq=vertcat(freql,freqh); sX2=vertcat(sX2l,sX2h); 
sY2=vertcat(sY2l,sY2h); cXY=vertcat(cXYl,cXYh);  sCR=vertcat(sCRl,sCRh);
figure, subplot(211), semilogx(freq,[smooth(dbm(FRFs)),smooth(dbm(FRFn))]);
title('FRF estimate'); ylabel('magnitude [dB]');
subplot(212), semilogx(freq,phs(FRFs));
xlabel('freq [Hz]'); ylabel('phase [deg]');

% STEP 3 - MOTOR
% parametric estimators comparison:
n=6;                        % model order of denominator polynomial
mh=4; ml=0;                 % model orders of numerator polynomial
r=(1:390);                  % frequency lines with information
relvar=1e-10;               % relative variation of costfunction (stop)
iter=4e2;                   % maximum number of iterations (stop)
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
figure
subplot(211), semilogx(freq,[dbm(FRFs),dbm(FRFls),dbm(FRFnls),...
               dbm(FRFml),dbm(FRFbtls),dbm(FRFgtls)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on, legend('FRF','LSE','NLS','MLE','BTLS','GTLS')
subplot(212), semilogx(freq,[phs(FRFs),phs(FRFls),phs(FRFnls),phs(FRFml)...
              ,phs(FRFbtls),phs(FRFgtls)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
FRF4m=FRFs; FRF4mn = FRFn; SYSm = SYSbtls;

% STEP 4 - TABLE
% select and simplify the system
[zp,pp,kp]=zpkdata(SYSt,'v'); pp=sort(pp);
for k=1:length(zp)
    if real(zp(k)) > 0; zp(k)=-real(zp(k))+1i*imag(zp(k)); end
end
zp=zp(2:end); n=8; kp=1.18e+15;
pp(1)=-0; pp(2)=-0; pp=pp(1:n);
SYStc = zpk(zp,pp,kp);
FRFt=squeeze(freqresp(SYSt,freq*2*pi));
FRFtc=squeeze(freqresp(SYStc,freq*2*pi));
figure, subplot(211)
semilogx(freq,[dbm(FRF4t),dbm(FRFt),dbm(FRFtc),dbm(FRFn)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on, legend('FRFs','FRFplnt','FRFcln','FRFn')
subplot(212),semilogx(freq,[phs(FRF4t),phs(FRFt),phs(FRFtc)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
disp('plant poles & zeros'); disp(pole(SYStc)); disp(zero(SYStc));
ppt = pp;

% STEP 4 - MOTOR
% select and simplify the system
[zp,pp,kp]=zpkdata(SYSm,'v'); pp=sort(pp);
for k=1:length(zp)
    if real(zp(k)) > 0, zp(k)=-real(zp(k))+1i*imag(zp(k)); end
end
n=6; pp(1)=-0; pp(2)=-0; pp=pp(1:n); 
SYSmc = zpk(zp,pp,kp);
FRFm=squeeze(freqresp(SYSm,freq*2*pi));
FRFmc=squeeze(freqresp(SYSmc,freq*2*pi));
figure, subplot(211)
semilogx(freq,[dbm(FRF4m),dbm(FRFm),dbm(FRFmc),dbm(FRF4mn)],'LineWidth',2)
ylabel('Amplitude [dB]'),grid on,legend('FRFs','FRFplnt','FRFcln','FRFn')
subplot(212),semilogx(freq,[phs(FRF4m),phs(FRFm),phs(FRFmc)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
disp('plant poles & zeros'); disp(pole(SYSmc)); disp(zero(SYSmc));
ppm = pp;

% STEP 5 - COMBINE: N4
% combine the design to SIMO and display
disp('plant poles'); disp(pole(SYSmc)); disp(pole(SYStc));
disp('plant zeros'); disp(zero(SYSmc)); disp(zero(SYStc));
figure, pzmap(SYSmc ,SYStc);
[zm,pm,km]=zpkdata(SYSmc,'v'); sys1=SYSmc;
[zl,pl,kl]=zpkdata(SYStc,'v'); sys2=zpk(zl,[pm;pl(5:6)],kl);
sys8 = ss([sys1;sys2]);
%[sysb, g] = balreal(sys8);
%sys8 = modred(sysb,[10,11,12,13,14],'del');
FRFm8=squeeze(freqresp(sys8(1),freq*2*pi));
FRFt8=squeeze(freqresp(sys8(2),freq*2*pi));
figure;
subplot(221); semilogx(freq,[dbm(FRF4m),dbm(FRFmc),dbm(FRFm8)],'LineWidth',2)
subplot(223); semilogx(freq,[phsw(FRF4m),phsw(FRFmc),phsw(FRFm8)],'LineWidth',2)
subplot(222); semilogx(freq,[dbm(FRF4t),dbm(FRFtc),dbm(FRFt8)],'LineWidth',2)
subplot(224); semilogx(freq,[phs(FRF4t),phs(FRFtc),phs(FRFt8)],'LineWidth',2)

% STEP 6 - FIGURE
figure;
subplot(221);A1=gca;P1=semilogx(freq,[dbm(FRF4m),dbm(FRF4mn),dbm(FRFm8)]);
subplot(223);A2=gca;P2=semilogx(freq,phsw(FRF4m),freq,phsw(FRFm8),'r');
subplot(222);A3=gca;P3=semilogx(freq,[dbm(FRF4t),dbm(FRF4tn),dbm(FRFt8)]);
subplot(224);A4=gca;P4=semilogx(freq,phs(FRF4t),freq,phs(FRFt8),'r');
S.Haxis=[A1,A2,A3,A4]; S.Hplot=[P1;P2;P3;P4];
S.Xlabel = {'','frequency [Hz]','','frequency [Hz]'};
S.Ylabel = {'magnitude [dB]','phase [deg]','',''};
S.Tlabel = {'motor-side','','table-side',''};
S.Llabel = {'measurement','noise model','plant model'};
S.Lpos = [0.43,0.42,0.2,0.2]; %lbwh
%S.Xlim = {0,sim_end;0,sim_end;0,sim_end;0,sim_end};
%S.Ylim = {0,1.01;0,max(xp(:,2))+1;-1e-4,1e-4;-1,1};
%S.Xtik = {0,0.5/df,1/df;10,1e2,1e3;0,tnrofs/fs/2,tnrofs/fs;10,1e2,1e3};
S.Oname = 'plot\modeln8';
%figpublish(S);

% FINAL
% save model and parameter estimations:
P=load('data\PsysData.mat');
P.sys40n8 = sys8; 
P.sys40n8t = SYStc; P.sys40n8m = SYSmc;
P.FRFs40h=[FRF4t,FRF4m]; P.FRFn40h=[FRF4tn,FRF4mn]; P.freq40h = freq;
save('data\PsysData.mat','-struct','P');