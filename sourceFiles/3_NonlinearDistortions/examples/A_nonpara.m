%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 FRF DATA                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_sf;
% Calculate FRF's
for d=1:nrofd
    [X(:,d),Y(:,d),freq,sX2(:,d),sY2(:,d),cXY(:,d),FRFms(:,d),sCRy(:,d)]...
                   = time2frfms(x(:,d),y(:,d),fs,nl,nh,nrofs);
    [Xn(:,d),Yn(:,d),sX2n(:,d),sY2n(:,d),cXYn(:,d),FRFmn(:,d)]...
                   = time2frfmn(x(:,d),y(:,d),fs,nl,nh,nrofs);
    snr_y(:,d)=20*log10(abs(Y(:,d))./sqrt(sY2(:,d)));
               
    [~,V(:,d),~,~,sV2(:,d),cXV(:,d),FRFvs(:,d),sCRv(:,d)]...
                   = time2frfvs(x(:,d),v(:,d),fs,nl,nh,nrofs);
    [~,Vn,~,sV2n(:,2),cXVn(:,d),FRFvn(:,d)]...
                   = time2frfmn(x(:,d),v(:,d),fs,nl,nh,nrofs);
    snr_v(:,k)=20*log10(abs(V(:,d))./sqrt(sV2(:,d)));
end
X=mean(X,2); Y=mean(Y,2); V=mean(V,2);
Xn=mean(Xn,2); Yn=mean(Yn,2); Vn=mean(Vn,2);
sX2=mean(sX2,2); sX2n=mean(sX2n,2);
sY2=mean(sY2,2); sY2n=mean(sY2n,2);
sV2=mean(sV2,2); sV2n=mean(sV2n,2);
cXY=mean(cXY,2); cXV=mean(cXV,2);
cXYn=mean(cXYn,2); cXVn=mean(cXVn,2);
sCRy=mean(sCRy,2); sCRv=mean(sCRv,2);
FRFms=mean(FRFms,2); FRFmn=mean(FRFmn,2); 
FRFvs=mean(FRFvs,2); FRFvn=mean(FRFvn,2);
snr_y=mean(snr_y,2); snr_v=mean(snr_v,2);
save(name,'X','Y','V','Xn','Yn','Vn','sX2','sX2n','sY2','sY2n',...
    'sV2','sV2n','cXY','cXYn','cXV','cXVn','sCRy','sCRv',...
    'FRFms','FRFmn','FRFvs','FRFvn','freq','snr_y','snr_v');

% Plot FRF results
%P_F1frfm;
%P_F2frfv;

% Plot SNR
P_F3snr;

%save('SNR','snr_y')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%