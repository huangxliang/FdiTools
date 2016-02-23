function [Xs,Ys,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,y,fs,fl,fh,nrofs)
%TIME2FRF_ML - maximum-likelihood estimation of FRF.
% [Xs,Ys,FRFs,FRFn,freq,sX2,sY2,cXY,sCR] = time2frf_ml(x,y,fs,fl,fh,nrofs)
% x, y      : data vectors measured using periodic broad-band excitation
% fs        : sampling frequency
% fl,fh     : lowest & highest frequency of excitated frequency band
% nrofs     : number of sample points per period
% FRFm,FRFn : measurement & noise frequency response functions
% sX2,sY2   : variance of real & imaginary parts of X,Y noise
% cXY       : covariance between real & imaginary parts of X,Y noise 
% sCR       : cramer-rao variance on measurement FRF
% Algorithm : non-parametric frequency domain identification (ch2)
%             System Identification: A Frequency Domain Approach,
%             R.Pintelon, J.Schoukens (2012), IEEE Press-Wiley.
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014

% PAR
x=x(:); y=y(:);                                   % input/output vectoring
fr = fs/nrofs;                                    % frequency resolution
nl = ceil(fl/fr); nh = floor(fh/fr);              % low & high freq points
freq = double((nl:1:nh)'/(nrofs/fs));             % freq lines
nr=length(freq);                                  % number of freq lines
nrofp = double(floor(length(x)/nrofs));           % number of periods

% FFT
INPs=[]; OUTs=[];
for i=1:nrofp
	INPs=[INPs fft(x(1+(i-1)*nrofs:i*nrofs))];    % fft of 1x period
	OUTs=[OUTs fft(y(1+(i-1)*nrofs:i*nrofs))];
end
INPs=INPs(nl+1:1:nh+1,:);                         % fft dc-term removal
OUTs=OUTs(nl+1:1:nh+1,:);

% VAR
sX2=((std(INPs'))'.^2)/2/nrofp;                   % measurement variances
sY2=((std(OUTs'))'.^2)/2/nrofp;
for (i=1:nr)
    Q=cov(INPs(i,:),OUTs(i,:));                   % measurement covariance
    cXY(i,1) = Q(1,2)/2/nrofp;
end

% FRF
Xs=mean(INPs,2); Ys=mean(OUTs,2);                 % frequency averaging
FRFs=Ys./Xs;                                      % signal frf calc
sCR=2*abs(FRFs).*(sX2./(abs(Xs)).^2 ...           % cramer-rao lowerbound
    +sY2./(abs(Ys)).^2 ...
    -2*real(cXY./(conj(Xs).*Ys)));

% NOISE
OUT2=[]; m=1;                             
for i=1:nrofp/2
 OUT2=[OUT2 fft(y(1+(i-1)*nrofs*2:i*nrofs*2))];   % fft of 2x period
end
OUT2=OUT2(2*nl:1:2*nh+1,:);                       % fft dc-term removal
    
for k=1:nr*2
    if mod(k,2)~=0                                % uneven freq lines
        OUTn(m,:)=OUT2(k,:);                      % noise data
        m=m+1;
    end
end
Yn=mean(OUTn,2);                                  % frequency averaging
FRFn = Yn./Xs;                                    % noise frf calc
end
