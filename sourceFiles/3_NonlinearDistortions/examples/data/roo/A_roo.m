%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              MULTISINE DATA                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Test parameters
load('P_roo.mat')
name='data\roo\D_roo';
fs=100e3;
nrofs=p*fs; df=fs/nrofs;
nl=ceil(fmin/df); nh=floor(fmax/df);
dt=1/fs; t=(0:dt:p-dt)';

% Test values
data(1)=load('T_roo1.mat'); data(2)=load('T_roo2.mat');
data(3)=load('T_roo3.mat'); data(4)=load('T_roo4.mat');
data(5)=load('T_roo5.mat'); nrofd=5;

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