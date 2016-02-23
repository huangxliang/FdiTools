%
% Calculation of full orthogonal random phase multisine excitations for
% measuring the frequency response matrix of a system with 3 inputs
%
% Notes: 
%
%       1. the same harmonics are excited in all inputs
%       2. the rms-values of the inputs may differ
%       3. the shape of the amplitude spectrum of each input can be chosen
%
% Rik Pintelon, October 2009
% 

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition spectral content multisine excitations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu = 3;                                     % number of inputs
N = 256*2;                                  % number of time domain samples 
ExcitedHarm = [N/8:1:N/4];                  % excited harmonics
nh = length(ExcitedHarm);

% definition amplitude spectra inputs
AmplitudeExcitedHarm = ones(nu, nh);
AmplitudeExcitedHarm(:, floor(nh/2):end) = 0.1; % relative amplitudes

% defition rms values inputs
RmsValues = [1; 10; 100];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of the nu x nu orthogonal multisine excitation signals %
% each time the function is called a new random phase is generated   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TheSignal = Calc_MIMO_Multisine(ExcitedHarm, N, AmplitudeExcitedHarm, RmsValues);

% rms values of the inputs
rms(TheSignal, 3)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the spectra of the orthogonal multisines %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = fft(TheSignal, [], 3)/sqrt(N);
rms(U, 3)

figure(1)
fig_num = 0;
for ii = 1:nu
    for jj = 1:nu
        fig_num = fig_num+1;
        subplot(nu, nu, fig_num);
        plot(db(squeeze(U(ii,jj,1:N/2))), '+')
        Old_axis = axis;
        axis([0, N/2, Old_axis(3:4)])
        zoom on
    end % jj column index
end % ii row index
subplot(nu, nu, 2)
title('DFT spectra')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether the phase matrix is unitary %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% at each excited harmonic the phase matrix of the nu x nu DFT spectra
% should be unitary

rel_error = ones(nh, 1);
Tphase_all = zeros(nu, nu, nh);
Tampl_all = zeros(nu, nu, nh);
for kk = 1:nh
    TheIndex = ExcitedHarm(kk)+1;
    T = squeeze(U(:,:,TheIndex));
    Tampl = abs(T)*sqrt(nu);
    Tphase = T./Tampl;
    rel_error(kk) = norm(Tphase' - inv(Tphase))/eps;
    Tphase_all(:,:,kk) = Tphase;
    Tampl_all(:,:,kk) = Tampl;
end % kk harmonic number

figure(2)
plot(ExcitedHarm, rel_error)
xlabel('harmonic number')
ylabel('relative error (in eps)')
title('check unitary phase matrix orthogonal multisines')

% show random phase distribution over the excited harmonics
figure(3)
figNum = 0;
for ii = 1:nu
    for jj = 1:nu
        figNum = figNum+1;
        subplot(nu, nu, figNum);
        plot(squeeze(Tphase_all(ii, jj, :))*sqrt(nu),'b+')
    end % jj column index
end % ii row index
subplot(nu, nu, 2);
title('Phase matrix DFT spectra')