%
% Single-input, single-output (SISO) system operating in open loop and excited by an random input
% Known input, noisy output case (generalized output error)
% Combining different experiments
%
% Rik Pintelon, 23 September 2011
%

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition discrete-time plant and noise models %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the discrete-time plant and noise models 
load PlantNoiseModel

N = 51200;                              % total number of time domain samples
N1 = 15000;                             % number of time domain samples in experiment no. 1
N2 = 25000;                             % number of time domain samples in experiment no. 2
N3 = N-(N1+N2);                         % number of time domain samples in experiment no. 3
fs = 19531;                             % sampling frequency
fstart = 1.9073e+001;                   % start frequency in Hz of the frequency band of interest
fstop = 5.8555e+003;                    % stop frequency in Hz of the frequency band of interest


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation known input - noisy output data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% experiment no. 1
u1 = randn(1, N1);                      % known input
y0_1 = filter(B, A, u1);                % noiseless output
e1 = randn(1, N1);                      % driving white noise source
v1 = filter(C, D, e1);                  % filtered white noise disturbance 
y1 = y0_1 + v1;                         % noisy output

% experiment no. 2
u2 = randn(1, N2);                      % known input
y0_2 = filter(B, A, u2);                % noiseless output
e2 = randn(1, N2);                      % driving white noise source
v2 = filter(C, D, e2);                  % filtered white noise disturbance
y2 = y0_2 + v2;                         % noisy output

% experiment no. 3
u3 = randn(1, N3);                      % known input
y0_3 = filter(B, A, u3);                % noiseless output
e3 = randn(1, N3);                      % driving white noise source
v3 = filter(C, D, e3);                  % filtered white noise disturbance
y3 = y0_3 + v3;                         % noisy output


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling the ArbLocalPolyAnal function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data 
data.u = [u1, u2, u3];                  % concatenation of the 3 experiments: row index is the input number; the column index the time instant (in samples) 
data.y = [y1, y2, y3];                  % concatenation of the 3 experiments: row index is the output number; the column index the time instant (in samples) 
data.Ts = 1/fs;                         % sampling period 
data.N = [N1, N2, N3];                  % number of samples in each experiment (in the same order as the concatenation)

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
title('Comparison estim. G and true G_0  FRF')
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

