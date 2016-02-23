function [hAx] = syssim(sys,t)
%SYSSIM typical input system simulation.
%
% sys   : simulated system
% author: Thomas Beauduin, KULeuven, 2014
%%%%%
% data extraction:
%sys=sys/freqresp(sys,0);

% UNFINISHED
% change impulse to zero-sin-zero
% add time analysis data on fig
% ex: settling time, overshoot, ...
% input signals:
imp = vertcat(1,zeros(size(t(2:end))));
stp = [zeros(round(0.1*size(t,1)),1);ones(round(0.9*size(t,1)),1)];
rmp = (0:(1/length(t)):1-(1/length(t)));
exp = t.^2./2;

figure
subplot(221),lsim(sys,imp,t,0); title('impuls response')
subplot(222),lsim(sys,stp,t,0); title('step response')
subplot(223),lsim(sys,rmp,t,0); title('ramp response')
subplot(224),lsim(sys,exp,t,0); title('exp response')
hAx = gca;

end

