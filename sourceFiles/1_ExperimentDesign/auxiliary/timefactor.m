function tf = timefactor(x,fmin,fmax,p)
%TIMEFACTOR - Time Factor (TF) of excitation signal.
%
% x         : time domain signal
% fmin,fmax : sample frequency
% Author    : Thomas Beauduin, KULeuven, 2015

% UNFINNISED

u = u(:);
n = length(u);
U = fft(u);
df=1/p;
nl=ceil(fmin/df); nh=floor(fmax/df);
N = nh-nl+1;

% Time factor
U_rms = sum(abs(U(2:N)).^2)/N;
tfv=0.5*cr^2*U_rms./(abs(U(2:N)).^2);
tf = max(tfv)

Armse = norm(U(nl+1:nh+1))/sqrt(N);
Amin = min(abs(U(nl+1:nh+1)));
tf2 = (0.5*((u_peak/Amin)^2)/2/N)*n*n
tf3 = 0.5*cr*cr*(Armse/Amin)^2

end
