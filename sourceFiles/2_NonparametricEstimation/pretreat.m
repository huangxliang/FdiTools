function [y,time] = pretreat(x,nrofs,fs,nroft)
%PRETREAT - pretreatment of data vectors for estimation.
% [y,time] = pretreat(x,nrofs,fs)
% x     : untreated periodic measurement vector
% nrofs : number of samples per signal period
% nroft : number of transient periods
% y     : transient, offset and trend removed data
% method: offset & trend removal from data
% Author: Thomas Beauduin, KULeuven, 2014

% 1. Transient Removal
x = x(nroft*nrofs+1:end);
nrofp = ceil(length(x)/nrofs);

% 2. Offset Removal
P = reshape(x,nrofs,nrofp);
Pa = mean(P,1);
for k=1:nrofp
    for l=1:nrofs
        P(l,k)=P(l,k)-Pa(k);
    end
end
y = P(:);
time = (0:(length(y)-1))'/fs;
end

% 3. Trend Removal:
% Peirlinckx, L., P. Guillaume, and R. Pintelon (1996). 
% Accurate and fast estimation of the Fourier coefficients of 
% periodic signals disturbed by a trend. 
% IEEE Trans. Instrum. and Meas., vol. 45, no.1, pp. 5-11.

% McCormack, A. S., J. O. Flower, and K. R. Godfrey (1994a). 
% The suppression of drift and transient effects for 
% frequecy-domain identification. 
% IEEE Trans. Instrum. and Meas., vol. 43, no. 2, pp.232-237.