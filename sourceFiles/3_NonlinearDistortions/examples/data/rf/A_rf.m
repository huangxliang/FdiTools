%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              MULTISINE DATA                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Test parameters
load('P_rf.mat');
name='data\rf\D_rf';
fs=100e3;
nrofs=p*fs; df=fs/nrofs;
nl=ceil(fmin/df); nh=floor(fmax/df);
dt=1/fs; t=(0:dt:p-dt)';

% Test values
data(1)=load('T_rfa1.mat'); data(2)=load('T_rfa2.mat');
data(3)=load('T_rfa3.mat'); data(4)=load('T_rfa4.mat');
data(5)=load('T_rfa5.mat');
data(6)=load('T_rfb1.mat'); data(7)=load('T_rfb2.mat');
data(8)=load('T_rfb3.mat'); data(9)=load('T_rfb4.mat');
data(10)=load('T_rfb5.mat');
data(11)=load('T_rfc1.mat'); data(12)=load('T_rfc2.mat');
data(13)=load('T_rfc3.mat'); data(14)=load('T_rfc4.mat');
data(15)=load('T_rfc5.mat');
data(16)=load('T_rfd1.mat'); data(17)=load('T_rfd2.mat');
data(18)=load('T_rfd3.mat'); data(19)=load('T_rfd4.mat');
data(20)=load('T_rfd5.mat');
nrofd=20;
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
    xi(d,:,:)=S'; yi(d,:,:)=A';
    x(:,d)=S(:); y(:,d)=A(:); v(:,d)=V(:);
end
clear V;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%