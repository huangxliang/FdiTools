%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           VWX MULTISINE DATA                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Test parameters
load('P_vwx3.mat');
name='data\vwx\D_vwx';
fs=100e3;
nrofs=p*fs; df=fs/nrofs;
nl=ceil(fmin/df); nh=floor(fmax/df);
dt=1/fs; t=(0:dt:p-dt)';

% Test values
data(1)=load('T_vwx3.mat');
nrofd=1;
for d=1:nrofd
    m=data(d);
    setp=m.UntitledSetPosition.Data;
    actp=m.UntitledActPosition.Data;
    vibp=m.UntitledVibrometer.Data;
    nrofp=length(setp)/nrofs;
 
    % Remove transient period
    S=reshape(setp,nrofs,nrofp); Sa=mean(S,1);
    A=reshape(actp,nrofs,nrofp); Aa=mean(A,1);
    V=reshape(vibp,nrofs,nrofp); Va=mean(V,1);
    S=S(:,2:end); A=A(:,2:end); V=V(:,2:end);
    nrofp=nrofp-1;

    % Remove steady-state offset
    for k=1:nrofp
        for l=1:nrofs
            S(l,k)=S(l,k)-Sa(k);
            A(l,k)=A(l,k)-Aa(k);
            V(l,k)=V(l,k)-Va(k);
        end
    end
    x(:,d)=S(:); y(:,d)=A(:); v(:,d)=V(:);
end
clear V;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%