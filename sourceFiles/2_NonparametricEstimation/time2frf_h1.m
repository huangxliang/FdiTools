function [FRF,freq,coh,X2,Y2] = time2frf_h1(x,y,fs,fl,fh,df,window)
%TIME2FRF_H1 - H1 estimation of FRF using a Hanning window.
%
% x, y     : unsynchonised input/output data vector
% fs       : measurement sampling frequency
% fl, fh   : lowest and highest frequencies of excitation
% nrofs    : number of sample points in one fft window
% window   : if 1 : hanning window, else no window
% FRF      : FRF-matrix [H11 H12 H13 .. H1m H21  ... Hl1 Hl2 ... Hlm]
% freq     : spectral lines of measured frf
% coh      : multiple coherence of frf
% X/Y      : output/input spectrum
% Author   : Thomas Beauduin, KULeuven, 2014
%
window = window(:);
[~,m] = size(x); 
[nroft,l] = size(y);

% Window Calculation
nrofs = fs/df;
nrofp = floor(nroft/nrofs);
if (window==1)
  win = hanning(nrofs);
else
  win = ones(nrofs,1);
end
winX = kron(win(:),ones(1,m));
winY = kron(win(:),ones(1,l));

% Windowed FFT
X = zeros(fh-fl+1,m*nrofp);
Y = zeros(fh-fl+1,m*nrofp);
for k=1:nrofp
  Xk = fft(x(1+(k-1)*nrofs:k*nrofs,:).*winX);
  Yk = fft(y(1+(k-1)*nrofs:k*nrofs,:).*winY);
  X(:,1+(k-1)*m:k*m) = Xk(fl+1:fh+1,:);
  Y(:,1+(k-1)*l:k*l) = Yk(fl+1:fh+1,:);
end

coh = zeros(fh-fl+1,l);
FRF = zeros(fh-fl+1,m*l);
freq= fs*(0:1:nrofs/2-1)'/nrofs;
freq= freq(fl/df+1:fh/df+1);

Xk = zeros(m,nrofp); 
Yk = zeros(l,nrofp);
Y2 = zeros(l*length(freq),l);
X2 = zeros(m*length(freq),m);

% FRF and COH
for k=1:fh-fl+1
  for i=1:nrofp
    Xk(:,i) = X(k,1+(i-1)*m:i*m).';
    Yk(:,i) = Y(k,1+(i-1)*l:i*l).';
  end
  XkXk = Xk*Xk';
  YkYk = Yk*Yk';
  
  H1kt = ((Yk*Xk')*inv(XkXk)).';
  H1k = H1kt.';
  FRF(k,:) = H1kt(:).';  
  for i=1:l
    coh(k,i) = real((H1k(i,:)*XkXk*H1k(i,:)')/YkYk(i,i));
  end
  X2(1+m*(k-1):m*k,:) = XkXk;
  Y2(1+l*(k-1):l*k,:) = YkYk;

end
end