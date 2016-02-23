function [y] = fullsignal(x,N)
%FULLSIGNAL - create full periodic signal.
%
% x         : single period of excitation signal 
% N         : total number of excitation points
% y         : repeated periodic excitation signal
% Author    : Thomas Beauduin, KULeuven
%             PMA division, February 2014
% see also MSIN, SSIN, PBRS

nrofs=length(x);            % number of samples
nrofp=floor(N/x);           % number of periods

y=zeros(nrofs*nrofp,1);
for k=0:nrofp-1
y(k*nrofs+1:(k+1)*nrofs,1)=x;
end

end

