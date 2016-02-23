%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 RAMP DATA                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
load('P_rp10.mat');
fs=100e3;
nrofs=round(p*fs);
t_rp10=(0:1/fs:(nrofs-1)/fs)';

% Configure data
m=load('T_rp10.mat');
set=m.UntitledSetPosition.Data;
act=m.UntitledActPosition.Data;
nrofp=round(length(set)/nrofs);
S=reshape(set,nrofs,nrofp);
A=reshape(act,nrofs,nrofp);
x0=mean(S,2); y0=mean(A,2);
offset=x0(100);
for l=1:length(x0);
    x_rp10(l,1)=x0(l)-offset;
    y_rp10(l,1)=y0(l)-offset;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%