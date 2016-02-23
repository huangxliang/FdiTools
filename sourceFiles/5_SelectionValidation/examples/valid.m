%
% create functions for the different tests
% 1 test = 1 function
% ex.: bundle cost_function calculation together
%
%% Cost function and interval
cost_nls=mlfdi_res(Bn_nls,An_nls,freq,X,Y,sX2,sY2,cXY,cORd,fs)
cost_ml=mlfdi_res(Bn_ml,An_ml,freq,X,Y,sX2,sY2,cXY,cORd,fs)
Vnoise =((nrofp-1)/(nrofp-2))*(length(freq)-n);
cost_intval=[Vnoise-2*Vnoise^(1/2) Vnoise+2*Vnoise^(1/2)]

%% Measured FRF - Modeled FRF
N_alfa=10.5966;
A1=smooth((abs(FRFms-FRFnls).^2)./1);
A2=smooth((abs(FRFms-FRFml).^2)./1);
B=N_alfa*sCRy/2;
figure
loglog(freq,[A1,A2,B])
title('Measured FRF - Modeled FRF vs. Cramer-Rao lower bound')
legend('NLSE','MLE','sCR')
ylabel('Error squared'), xlabel('frequency [Hz]')
ylim([10^-5 10^0])

%% Resid
[Lags,Corr(:,1),CB,Fraction(:,1)]=resid(freq,FRFms,FRFnls,sCRy,nrofp);
[Lags,Corr(:,2),CB,Fraction(:,2)]=resid(freq,FRFms,FRFml,sCRy,nrofp);
figure
plot(Lags,Corr,'*',Lags,CB(1,:),'k',Lags,CB(2,:),'k--')
title(strcat('Frac_{50%} :',num2str(Fraction(:,1)),...
             'Frac_{95%} :',num2str(Fraction(:,2))));
xlabel('{\itk}'), ylabel('{\itR} {\itk}')
legend('NLS','ML')

%% Simulation JUMP
A_rp10;
SYSnls=SYSnls/freqresp(SYSnls,0); SYSml=SYSml/freqresp(SYSml,0);
y_nls=lsim(SYSnls,x_rp10,t_rp10); y_ml=lsim(SYSml,x_rp10,t_rp10);
figure
plot(t_rp10,[x_rp10,y_rp10,y_nls,y_ml]);
legend('IN','OUT','NLS','ML')
title('V JUMP')

%% Simulation SCAN
A_rp2;
SYSnls=SYSnls/freqresp(SYSnls,0); SYSml=SYSml/freqresp(SYSml,0);
y_nls=lsim(SYSnls,x_rp2,t_rp2); y_ml=lsim(SYSml,x_rp2,t_rp2);
figure
plot(t_rp2,[x_rp2,y_rp2,y_nls,y_ml]);
legend('IN','OUT','NLS','ML')
title('V SCAN')

%save('S_fdi.mat','SYSml')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%