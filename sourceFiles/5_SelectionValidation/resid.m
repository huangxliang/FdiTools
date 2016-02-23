function [Lags,Corr,CB,Fraction] = resid(freq,FRFm,FRF,sCR,nrofp)
% RESID - Identification Residuals for validation
%
% freq      : freq vector of measurement
% FRFm      : FRF measurement
% FRF       : FRF of fitted model
% sCR       : Cramer-Rao underbound of measurement
% norfp     : nr of periods in measurment
% Lags      : lags for auto_corr ploting
% Corr      : Auto correlation
% CB        : Confindence Bounds
% Fraction  : fraction under confidence bounde
% Author    : Thomas Beauduin, KULeuven, 2014
%
F=length(freq);
Res=zeros(1,F);
Fraction = zeros(1,2);
Lags = (-F+1:1:F-1);
% Residuals
for kk = 1:F    
    Res(kk) = (FRFm(kk)-FRF(kk))./sCR(kk).^0.5;    
end
% auto-correlation
Auto_Corr = xcorr(Res,'unbiased')./F;
% std auto-correlation
scale0 = (nrofp-2)/(nrofp-1);
scale = (nrofp-5/3)/(nrofp-11/12);
Conf_scale0=scale0*((nrofp-1)^(3/2)/(nrofp-2)/(nrofp-3)^(1/2));
Conf_scale=scale*(nrofp-1)/(nrofp-2);
p50 = sqrt(-log(1-0.5));
p95 = sqrt(-log(1-0.95));
ConfScale = Conf_scale*ones(size(Lags));
ConfScale(F) = Conf_scale0;
Conf_Bound = repmat(ConfScale./(F-abs(Lags)).^0.5,[2,1]);
Conf_Bound(1,:) = p50*Conf_Bound(1,:);
Conf_Bound(2,:) = p95*Conf_Bound(2,:);
CB=Conf_Bound;

% fractions outside p% confidence bounds
Select = (1:1:2*F-1);
Select(F) = [];

Fraction(1) = length(find(abs(Auto_Corr(Select)) - ...
                          Conf_Bound(1,Select) > 0)) / (2*F-2)*100;
Fraction(2) = length(find(abs(Auto_Corr(Select)) - ...
                          Conf_Bound(2,Select) > 0)) / (2*F-2)*100;

Corr=squeeze(abs(Auto_Corr(:)));
% PLOT AUTO-CORRELATION
% figure 11-8
end

