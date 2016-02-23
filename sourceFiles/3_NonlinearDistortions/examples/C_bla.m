%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        BEST LINEAR APPROXIMATION                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_rf;

% Calculate BLA
for k=1:nrofd
    [X(:,k),Y(:,k),freq,sX2(:,k),sY2(:,k),cXY(:,k),FRFs(:,k),sCRy(:,k)]...
                     = time2frfms(x(:,k),y(:,k),fs,nl,nh,nrofs);
    [Xn,Yn,~,sY2n,cXYn,FRFn(:,k)]...
                     = time2frfmn(x(:,k),y(:,k),fs,nl,nh,nrofs);
end
FRFs=FRFs.*0.86; Y=Y.*0.86;
X=mean(X(10:500,:),2); Y=mean(Y(10:500,:),2); 
sX2=mean(sX2(10:500,:),2); sY2=mean(sY2(10:500,:),2);
cXY=mean(cXY(10:500,:),2); sCRy=mean(sCRy(10:500,:),2); 
FRFms=mean(FRFs(10:500,:),2); FRFmn=mean(FRFn(10:500,:),2);
freq=freq(10:500);
save(name,'X','Y','freq','sX2','sY2','cXY','sCRy','FRFms','FRFmn');

% Plot BLA measurement
%P_F5bla1;

% Visualize BLA
nroff=length(freq);
FRF=zeros(nrofd,nrofp,nroff);
for M=1:nrofd
    for P=1:nrofp
        [~,~,~,~,~,~,FRF(M,P,:),~]...
                  = time2frfms(xi(M,P,:),yi(M,P,:),fs,10,nh,nrofs);
    end
end
Gmean = squeeze(mean(mean(FRF,2),1)).*0.86;
Gstdt = squeeze(std(mean(FRF, 2), 0, 1))/sqrt(nrofd);
Gstdn = (squeeze(mean(std(FRF, 0, 2).^2, 1))/(nrofd*nrofp)).^0.5;
Gstds = (nrofd*abs(Gstdt.^2 - Gstdn.^2)).^0.5;

% Plot BLA visualization
F_Lbla;
P_F5bla2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Control correctness
% figure
% semilogx(freq,dbm(Gmean),'k','LineWidth',4), hold on 
% semilogx(freq,dbm(FRFs),'g',freq,dbm(FRFn),'r',freq,dbm(Gstdn), 'b');
% legend('BLA','BLA','stdn','stdn')