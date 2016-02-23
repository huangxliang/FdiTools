%
% Single-input, single-output (SISO) system operating in open loop and excited by an random input
% noisy input, noisy output case (errors-in-variables) 
%
% Rik Pintelon, 23 September 2011
%

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition discrete-time plant, actuator and noise models %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the discrete-time plant and output noise models 
load PlantNoiseModel

% actuator characteristics
[Bact, Aact] = butter(6, 2*0.3);

% input noise model
[Cin, Din] = cheby1(8, 15, 2*0.3);
Cin = 0.1*Cin;

N = 51200;                              % number of time domain samples
fs = 19531;                             % sampling frequency
fstart = 1.9073e+001;                   % start frequency in Hz of the frequency band of interest
fstop = 5.8555e+003;                    % stop frequency in Hz of the frequency band of interest


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation known input - noisy output data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reference signal
r = randn(1, N);

% deterministic part
u0 = filter(Bact, Aact, r);             % noiseless input
y0 = filter(B, A, u0);                  % noiseless output

% input-output noise disturbances
eu = randn(1, N);                       % driving white noise source
vu = filter(Cin, Din, eu);              % input noise 
ey = randn(1, N);                       % driving white noise source
vy = filter(C, D, ey);                  % output noise 

% noisy input-output
u = u0 + vu;
y = y0 + vy;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling the ArbLocalPolyAnal function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data 
data.u = u;                             % row index is the input number; the column index the time instant (in samples) 
data.y = y;                             % row index is the output number; the column index the time instant (in samples)
data.r = r;                             % known reference signal
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
varVy = squeeze(CZ.n(1,1,:));           % estimate output noise variance
varVu = squeeze(CZ.n(2,2,:));           % estimate input noise variance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison estimates and true values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% true FRF
q = exp(-sqrt(-1)*2*pi*freq/fs);        % z^(-1) as a function of the frequency 
G0 = polyval(fliplr(B), q) ./ polyval(fliplr(A), q);
G0 = G0.';

% true input noise model and input variance
Hu_0 = polyval(fliplr(Cin), q) ./ polyval(fliplr(Din), q);
varVu_0 = abs(Hu_0.^2).';                 % true output noise variance

% true output noise model and output variance
Hy_0 = polyval(fliplr(C), q) ./ polyval(fliplr(D), q);
varVy_0 = abs(Hy_0.^2).';                 % true output noise variance

% plot the true plant and input-output noise models
figure(1)
plot(freq, db(G0), 'k', freq, db(Hy_0), 'r', freq, db(Hu_0), 'b')
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
title('True plant G_0, input H_{u_0} and output H_{y_0} noise models')
legend('G_0', 'H_{y_0}', 'H_{u_0}', 'Location', 'EastOutside');
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

% comparison true and estimated input-output noise variances
figure(3)
subplot(211)
plot(freq, db(varVu)/2, 'r', freq, db(varVu_0)/2, 'k');
xlabel('Frequency (Hz)')
ylabel('Input variance (dB)')
legend('estim. input var.', 'true input var.', 'Location', 'EastOutside');
subplot(212)
plot(freq, db(varVy)/2, 'r', freq, db(varVy_0)/2, 'k');
xlabel('Frequency (Hz)')
ylabel('Output variance (dB)')
legend('estim. output var.', 'true output var.', 'Location', 'EastOutside');
zoom on;
shg

