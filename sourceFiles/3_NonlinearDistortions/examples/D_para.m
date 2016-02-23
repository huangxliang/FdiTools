%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PARAMETRIC FDI                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_bla; close all;
n=3;
mh=1;
ml=0;
r=(1:300);

%%%%%
FRF_W = ones(size(FRFms)); iter=1000;
relvar=1e-10; GN = 0; cORd = 'c';

[Bn_nls,An_nls,~,~,Bls2,Als2] ...
= nllsfdi(FRFms(r),freq(r),FRF_W(r),n,mh,ml,iter,relvar,GN,cORd,fs);
SYSls = tf(Bls2,Als2);
FRFls=squeeze(freqresp(SYSls,freq*2*pi));
SYSnls = tf(Bn_nls,An_nls);

[Bn_ml,An_ml,Bls,Als,cost0,costls] ...
= mlfdi(X(r),Y(r),freq(r),sX2(r),sY2(r),cXY(r),n,mh,ml,iter,relvar,cORd,fs);
SYSml = tf(Bn_ml,An_ml);

% Alter system (remove NMP - continu model)
[znls,pnls,knls]=zpkdata(SYSnls); [zml,pml,kml] = zpkdata(SYSml);
znls=cell2mat(znls); zml=cell2mat(zml);
for k=1:mh
    if real(znls(k)) > 0
        znls(k)=-real(znls(k))+1i*imag(znls(k));
        zml(k)=-real(zml(k))+1i*imag(zml(k));
    end
end
SYSnls = zpk(znls,pnls,knls);
SYSml=zpk(zml,pml,kml);

%SYSnls=SYSnls/freqresp(SYSnls,0);
%SYSml=SYSml/freqresp(SYSml,0);

FRFnls=squeeze(freqresp(SYSnls,freq*2*pi));
FRFml=squeeze(freqresp(SYSml,freq*2*pi));

% Plot estimation result
figure
subplot(211)
semilogx(freq,[dbm(FRFms) dbm(FRFnls) dbm(FRFml)],'LineWidth',2)
ylabel('Amplitude [dB]'), grid on
legend('FRF','NLSE','MLE')
xlim([100,5000])
subplot(212)
semilogx(freq,[phs(FRFms) phs(FRFnls) phs(FRFml)],'LineWidth',2)
ylabel('Phase [deg]'),xlabel('frequency [Hz]'),grid on
xlim([100,5000])

% Plot pole-zero plot
figure
pzmap(SYSnls,SYSml)
SYSnlsd=c2d(SYSnls,1/(200e3));
SYSmld=c2d(SYSml,1/(200e3));
figure
pzmap(SYSnlsd,SYSmld)

% Alter system (remove NMP - discrete model)
[znls,pnls,knls]=zpkdata(SYSnlsd); [zml,pml,kml] = zpkdata(SYSmld);
znls=cell2mat(znls); zml=cell2mat(zml);
for k=1:mh
    if real(znls(k)) > 0
        znls(k)=-real(znls(k))+1i*imag(znls(k));
        zml(k)=-real(zml(k))+1i*imag(zml(k));
    end
end
SYSnlsd2 = zpk(znls,pnls,knls);
SYSmld2=zpk(zml,pml,kml);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%