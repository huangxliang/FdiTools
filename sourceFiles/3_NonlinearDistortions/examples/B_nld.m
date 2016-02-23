%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         NON-LINEAR DISTORTIONS                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_roo;

%Calculate FRF
for k=1:nrofd
    [X,Y(:,k),freq,sX2,sY2,cXY,FRFms,sCRy]...
                   = time2frfms(x(:,k),y(:,k),fs,nl,nh,nrofs);
    [Xn,Yn(:,k),~,sY2n,cXYn,FRFmn]...
                   = time2frfmn(x(:,k),y(:,k),fs,nl,nh,nrofs);
end
Y=mean(Y,2); Yn=mean(Yn,2); 

% FRF + odd distortions
freqA=[]; FRFA=[];
for i = 1:4:length(freq) 
    freqA = vertcat(freqA,freq(i));
    FRFA = vertcat(FRFA,Y(i)); 
end
% Even distortions
freqB=[]; FRFB=[];
for i = 2:2:length(freq)
    freqB = vertcat(freqB,freq(i));
    FRFB = vertcat(FRFB,Y(i)); 
end
% Odd distortions
freqC=[]; FRFC=[];
for i = 3:4:length(freq)
    freqC = vertcat(freqC,freq(i));
    FRFC = vertcat(FRFC,Y(i));
end
FRFC=FRFC.*1.5;
% save(name,'freqA','freqB','freqC','freq','FRFA','FRFB','FRFC','Yn');
% 
% A=load('D_roo.mat'); B=load('D_roo.mat');
% freqA=vertcat(A.freqA,B.freqC); freqB=vertcat(A.freqB,B.freqB);
% freqC=vertcat(A.freqC,B.freqA); freq=vertcat(A.freq,B.freq);
% FRFA=vertcat(A.FRFA,B.FRFC); FRFB=vertcat(A.FRFB,B.FRFB);
% FRFC=vertcat(A.FRFC,B.FRFA); Yn=vertcat(A.Yn,B.Yn);
% Plot Results
P_F4nld;
F_Lnld;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%