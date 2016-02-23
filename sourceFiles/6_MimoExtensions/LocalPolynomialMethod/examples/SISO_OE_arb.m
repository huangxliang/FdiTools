%
% Single-input, single-output (SISO) system operating in open loop and excited by an random input
% Known input, noisy output case (generalized output error) 
%
% Rik Pintelon, 22 September 2011
%

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition discrete-time plant and noise models %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the discrete-time plant and noise models 
load PlantNoiseModel

N = 51200;                              % number of time domain samples
fs = 19531;                             % sampling frequency
fstart = 1.9073e+001;                   % start frequency in Hz of the frequency band of interest
fstop = 5.8555e+003;                    % stop frequency in Hz of the frequency band of interest


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation known input - noisy output data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deterministic part
u0 = randn(1, N);                       % known input
y0 = filter(B, A, u0);                  % noiseless output

% noise disturbance
e = randn(1, N);                        % driving white noise source
v = filter(C, D, e);                    % filtered white noise disturbance 

% noisy output
y = y0 + v;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling the ArbLocalPolyAnal function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data 
data.u = u0;                            % row index is the input number; the column index the time instant (in samples) 
data.y = y;                             % row index is the output number; the column index the time instant (in samples) 
data.Ts = 1/fs;                         % sampling period

% method
method.dof = 6;                         % degrees of freedom of the variance estimate
method.order = 2;                       % order local polynomial approximation
method.startfreq = fstart;              % defines the start frequency of the analysis 
method.stopfreq = fstop;                % defines the stop frequency of the analysis

% local polynomial estimate FRF and its variance
[CZ, Z, freq, G, CvecG, dof, CL] = ArbLocalPolyAnal(data, method);

G = squeeze(G);                         % FRF estimate
varG = squeeze(CvecG);                  % variance FRF estimate
varV = squeeze(CZ.n(1,1,:));            % estimate output noise variance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison estimates and true values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% true FRF
q = exp(-sqrt(-1)*2*pi*freq/fs);        % z^(-1) as a function of the frequency 
G0 = polyval(fliplr(B), q) ./ polyval(fliplr(A), q);
G0 = G0.';

% true noise model and output variance
H0 = polyval(fliplr(C), q) ./ polyval(fliplr(D), q);
varV0 = abs(H0.^2).';                     % true output noise variance

% plot the true plant and noise models
figure(1)
plot(freq, db(G0), 'k', freq, db(H0), 'k--')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
title('True plant G_0 and noise H_0 models')
legend('G_0', 'H_0', 'Location', 'EastOutside');
zoom on;
shg

% comparison true and estimated FRF
figure(2)
plot(freq, db(G), 'r', freq, db(G0), 'k', freq, db(G-G0), 'k--', freq, db(varG)/2, 'r--');
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
title('Comparison estimated G and true G_0  FRF')
legend('G-estimate', 'G_0', '|G-G_0|', 'var(G)', 'Location', 'EastOutside');
zoom on;
shg

% comparison true and estimated output noise variance
figure(3)
plot(freq, db(varV)/2, 'r', freq, db(varV0)/2, 'k');
xlabel('Frequency (Hz)')
ylabel('Variance (dB)')
title('Comparison estimated and true noise variance')
legend('estim. output var.', 'true output var.', 'Location', 'EastOutside');
zoom on;
shg

