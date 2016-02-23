%
% Single-input, single-output (SISO) discrete-time Wiener-Hammerstein (WH) system
% excited by a random phase multisines
% The discrete-time Wiener-Hammerstein system consists of the cascade of 
%  
%   1. a second order system:           G1(z) 
%   2. a static nonlinear function:     x + alfa*x^2 + beta*x^3 
%   3. a first order system:            G2(z)  
%
% Noisy input, noisy output case (errors-in-variables problem) 
% Note: this examples uses the "DesignMultisineExcitations" toolbox
%
% Rik Pintelon, 4 October 2011
%

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition multisine excitation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExcitedHarm = [1401:1:2200];            % excited harmonics multisine excitation
fs = 4e9;                               % sampling frequency
N = 8000;                               % number of samples in one period
rms_r = 1.5;                            % rms value reference signal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disturbing input/output noise %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input noise
stdu = 0.2;      % input noise std
[c_in, d_in] = cheby1(4, 6, 2*[ExcitedHarm(1)/N, ExcitedHarm(end)/N]);

% output noise
stdy = 0.3;      % output noise std
[c_out, d_out] = cheby1(6, 10, 2*[ExcitedHarm(1)/N, ExcitedHarm(end)/N]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition discrete-time Wiener-Hammerstein system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first linear dynamic block
a1 = [1, -0.2, 0.9];
b1 = [1];

% static nonlinearity z = x + alpha*x^2 + beta*x^3
alpha = 0.1;
beta = 0.001;

% second linear dynamic block
a2 = [1, -0.5, 0.9];
b2 = [0, 1, 0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition actuator characteristics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b_act, a_act] = butter(6, 2*0.28);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transient response WH-system to random phase multisine excitations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 2;              % number of realisations random phase multisine
P = 2;              % number of consecutive periods of the transient response

% data structure
data = struct('u', [], 'y', [], 'r', [], 'N', [], 'Ts', [], 'ExcitedHarm', []);
data.u = zeros(1, 1, M, N*P);       % input signals
data.y = zeros(1, 1, M, N*P);       % ouput signals
data.r = zeros(1, 1, M, N);         % reference signals (one period of each realisation)
data.N = N;                         % samples in one period
data.Ts = 1/fs;                     % sampling period
data.ExcitedHarm = ExcitedHarm;     % excited harmonics multisine

% transient response to M different random phase realisations multisine excitation
for mm = 1:M
    
    % 1 period random phase multisine r(t) with rms value equal to r_rms
    r = CalcMultisine(ExcitedHarm, N);  % rms value = 1; equal harmonic amplitudes
    r = rms_r*r.';
    
    % transient response actuator over P periods
    u0 = filter(b_act, a_act, repmat(r, [1, P]));
    
    % transient response first LTI block
    x = filter(b1, a1, u0);
    
    % static nonlinearity
    z = x + alpha*x.^2 + beta*x.^3;
    
    % transient response second LTI block
    y0 = filter(b2, a2, z);
    
    % noisy input observation
    eu = randn(1, P*N);                             % white noise
    nu = stdu * filter(c_in, d_in, eu);             % filtered white input noise
    data.u(1,1,mm,:) = u0 + nu;
    
    % noisy output observation
    ey = randn(1, P*N);                             % white noise
    ny = stdy * filter(c_out, d_out, ey);           % filtered white output noise
    data.y(1,1,mm,:) = y0 + ny;
    
    % reference signal
    data.r(1,1,mm,:) = r;
      
end % for


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling the RobustLocalPolyAnal function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method.dof = 6;                         % degrees of freedom variance estimates

[CZ, Z, freq, G, CvecG, dof, CL] = RobustLocalPolyAnal(data, method);

% FRF and its variance
G = squeeze(G).';
varGn = squeeze(CvecG.n).';             % noise variance
varGNL = squeeze(CvecG.NL).';           % total variance (noise + NL distortions)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True best linear approximation (BLA)  = gamma * Gfirst * Gsecond %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = exp(-sqrt(-1)*2*pi*freq/fs);                        % powers of z^-1
G1 = polyval(fliplr(b1), q)./polyval(fliplr(a1), q);    % FRF first LTI block
G2 = polyval(fliplr(b2), q)./polyval(fliplr(a2), q);    % FRF second LTI block
G_BLA = G1.*G2;
% estimate gamma in least square sense
gamma = sum(real(G_BLA .* conj(G))) / sum(abs(G_BLA.^2));
G_BLA = gamma*G_BLA;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True input-output noise variances %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input noise variance
varU0 = abs(stdu*polyval(fliplr(c_in), q)./polyval(fliplr(d_in), q)).^2;
varU0 = varU0 / (M*P*N);    % averaging over P periods, M realisations + scaling N Fourier coefficients

% output noise variance
varY0 = abs(stdy*polyval(fliplr(c_out), q)./polyval(fliplr(d_out), q)).^2;
varY0 = varY0 / (M*P*N);    % averaging over P periods, M realisations + scaling N Fourier coefficients


%%%%%%%%%%%
% Figures %
%%%%%%%%%%%

freq_s = freq/1e9;      % scaled frequency

% comparison true and estimated input-output noise variances
% true input noise var = stdu^2 / (P*N): averaging over P periods + scaling by N for Fourier coefficients 
% true output noise var = stdy^2 / (P*N): averaging over P periods + scaling by N for Fourier coefficients 
% Note: due to the transient removal the estimated input-output noise variances are about 1 dB above the true value 
figure(1)
subplot(211)
plot(freq, db(squeeze(CZ.n(2,2,1,:)))/2, 'g', freq, db(varU0)/2, 'k');
xlabel('Frequency (Hz)')
ylabel('Input noise variance (dB)')
legend('estimate', 'true value', 'Location', 'EastOutside');
subplot(212)
plot(freq, db(squeeze(CZ.n(1,1,1,:)))/2, 'g', freq, db(varY0)/2, 'k');
xlabel('Frequency (Hz)')
ylabel('Output noise variance (dB)')
legend('estimate', 'true value', 'Location', 'EastOutside');

% estimated BLA, its noise and total variances
figure(2)
plot(freq_s, db(G), 'k', freq_s, db(varGn)/2, 'g', freq_s, db(varGNL)/2, 'r')
xlabel('Frequency (GHz)')
ylabel('{\itG}_{BLA} (dB)')
title('Estimated BLA and its variances')
legend('{\itG}_{BLA}', 'noise variance', 'total variance', 'Location', 'EastOutside');
zoom on
shg

% comparison true BLA with estimated BLA
figure(3)
plot(freq_s, db(G), 'k', freq_s, db(G_BLA), 'k--', ...
     freq_s, db(G-G_BLA), 'b', freq_s, db(varGNL)/2, 'r')
xlabel('Frequency (GHz)')
ylabel('Amplitude (dB)')
title('Comparison estimated and true BLA')
legend('estimated BLA', 'true BLA', 'difference', 'total variance', 'Location', 'EastOutside');
zoom on
shg
