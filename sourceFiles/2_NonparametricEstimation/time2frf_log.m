function [FRF,freq,coh] = time2frf_log(x,y,fs,fl,fh,df,window,nrofl)
%TIME2FRF_LOG - logarithmic estimation of frf.
%
% x, y     : input and output measured data vector
% fs       : measurement sampling frequency
% fl, fh   : lowest and highest frequencies of excitation
% df       : spectral line differences
% window   : 
% nrofl    : number of samples for window overlap
% FRF      : FRF-matrix [H11 H12 H13 .. H1m H21 ... Hlm]
% freq     : measured frequency lines
% coh      : multiple coherence
% Author   : Thomas Beauduin, KULeuven, 2015
%
window = window(:);
[~,m] = size(x); 
[nroft,l] = size(y);

% Windowed FFT
nrofw = length(window);
if nroft < nrofw
    x(nrofw)=0;  nroft=nrofw;
end
k = fix((nroft-nrofl)/(nrofw-nrofl))
index = (1:nrofw);
X = zeros(fh-fl+1,m*k); Y = zeros(fh-fl+1,l*k);
for i=1:k
   xw = kron(ones(1,m),window).*(x(index,:));
   yw = kron(ones(1,l),window).*(y(index,:));
   Xw = fft(xw); Yw = fft(yw);
   X(:,1+(i-1)*m:i*m) = Xw(fl+1:fh+1,:);
   Y(:,1+(i-1)*l:i*l) = Yw(fl+1:fh+1,:);
   index = index + (nrofw - nrofl);
end

% FRF Calculations
freq = (fl:fh)'*df;
FRF = zeros((fh-fl+1),m*l);
coh = zeros(fh-fl+1,l);
Xi = zeros(m,m);Yi=zeros(l,m);
for (fk = 1:(fh-fl+1))
  Hi = zeros(l,m);
  for (i=1:k-m+1)
     for (ii=1:m)
       Xi(:,ii) = X(fk,1+(ii+i-2)*m:(ii+i-1)*m).'; 
       Yi(:,ii) = Y(fk,1+(ii+i-2)*l:(ii+i-1)*l).'; 
     end
     Hi =  Hi + Yi*inv(Xi);
  end
  Hi = unwrap(angle(Hi/(k-m+1)));
  Hav_fk = zeros(l,m);
  for (i=1:k-m+1)
    for (ii=1:m)
       Xi(:,ii) = X(fk,1+(ii+i-2)*m:(ii+i-1)*m).'; 
       Yi(:,ii) = Y(fk,1+(ii+i-2)*l:(ii+i-1)*l).'; 
    end
    Hav_fk = Hav_fk + log(Yi*inv(Xi).*exp(-1i*Hi));
  end
  Hav_fk = exp(Hav_fk/(k-m+1)).*exp(1i*Hi); 
  Hav_fk = Hav_fk.'; Hav_fk = Hav_fk(:);
  FRF(fk,:) = Hav_fk.';
  XiXi = Xi*Xi'; YiYi = Yi*Yi';
  for i=1:l
    coh(fk,i) = real((Hav_fk(i,:)*XiXi*Hav_fk(i,:)')/YiYi(i,i));
  end
end
end