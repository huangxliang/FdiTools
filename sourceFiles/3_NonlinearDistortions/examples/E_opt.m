%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       OPTIMAL EXPERIMENT DESIGN                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D_para;

Fipp = cr_rao(X,Y,freq,Bn_nls,An_nls,sX2,sY2,cXY,n,mh,ml,cORd,fs);
for k=1:length(freq)
    Fippk=cr_rao(1,Y(k),freq(k),Bn_nls,An_nls,sX2(k),...
                  sY2(k),cXY(k),n,mh,ml,cORd,fs);
    vwx(k) = trace(inv(Fipp).*Fippk);
end

% Plot dispersion function
P_F6opt;

save('vwxm','vwx');
% sys = frd(vwx+1,freq.*(pi*2));
% b1 = fitmagfrd(sys,3);
% b1g = frd(b1,freq.*(2*pi));
% figure
% bodemag(sys,'r',b1g,'k:');
% 
% save('vwx','b1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%