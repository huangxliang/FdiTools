% read me file: use of the functions for designing random phase multisine excitations 
% Examples illustrating the use of the functions can be found in the folder "examples" 
%
%%	1. Design random phase multisines with random harmonic grid 
%
%       Two types of multisines are considered:
%
%           1. Odd multisines: excite the odd harmonics only
%           2. Full multisines: excite the even and odd harmonics
%
%       The harmonics can have the following frequency spacing:
%
%           1. Linear spacing
%           2. Logarithmic spacing (rounded to the nearest DFT line)
%
%       To detect the nonlinear distortions, 1 out of Nblock consecutive excited harmonics is
%       randomly eliminated. The harmonic grid satisfying the previous constraints is calculated by the function:
%
%           [ExcitedHarm, N, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, TypeMulti); 
%
%       Given the excited harmonics, the time signal is calculated using:
%
%           TheSignal = CalcMultisine(ExcitedHarm, N, AmplExcitedHarm);
%
%       Each time the function is called another random realisation of the phases is generated. 
%
%       Usage for single input systems:
%           1. Odd random phase multisines with random harmonic grid:
%              detection and classification in even and odd distortions.
%           2. Full random phase multisines with random harmonic grid:
%              detection of the total (even + odd) nonlinear distortions.
%           3. Measurement of the frequency response function, its noise 
%              variance and its total (noise + nonlinear distortions) variance. 
%
%       References:
%
%           Vanhoenacker, K., T. Dobrowiecki, and J. Schoukens (2001). Design of multisine excitations to characterize 
%           the nonlinear distortion during FRF-measurements, IEEE Trans. Instrum. Meas., vol. 50, no. 5, pp. 1097-1102. 
%
%           Pintelon R., G. Vandersteen, L. De Locht, Y. Rolain and J. Schoukens (2004). Experimental characterization of 
%           operational amplifiers: a system identification approach - Part I: theory and simulations, IEEE Trans. Instrum. Meas., 
%           vol. 53, no. 3, pp. 854-862.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA).
%
%%  2. Design of full random orthogonal multisines for multi-input systems
%
%       For each of the nu MIMO experiments, the nu random phase multisines excite the same harmonics. 
%       The nu sets of nu multisines are designed in such way that their DFT spectra are orthogonal:
%
%           1. At each excited harmonic the phase matrix of the nu x nu DFT spectrum is unitary. 
%           2. The rms values may differ over the inputs, but remain the same over the different 
%              MIMO experiments. 
%           3. The amplitude spectra may differ of the inputs, but remain the same over the different 
%              MIMO experiments.
%
%       Given the excited harmonics, the set of nu x nu full random orthogonal multisines are calculated 
%       by the function:
%
%           TheSignal = Calc_MIMO_Multisine(ExcitedHarm, N, AmplExcitedHarm, RmsValues); 
%
%       Each time the function is called another random realisation of the phases is generated. 
%
%
%       Usage for nu input systems:
%           1. Measurement of the frequency response matrix (FRM) and its
%           noise covariance. 
%           2. Measurement of the FRM, its noise covariance and its total noise + nonlinear distortions) covariance. 
%
%       References:
%
%           Dobrowiecki, T. P., J. Schoukens, and P. Guillaume (2006). Optimized excitation signals for MIMO frequency 
%           response function measurements, IEEE Trans. Instrum. and Meas., vol. 55, no. 6, pp. 2072–2079.
%
%           Dobrowiecki, T. P., and J. Schoukens (2007). Measuring a linear approximation to weakly nonlinear MIMO systems, 
%           Automatica, vol. 43, no. 10, pp. 1737–1751.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA).
%
%
%%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 4 October 2011 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%
