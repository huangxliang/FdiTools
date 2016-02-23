% Read me first file: use of the "FrequencyDomainToolbox":
% 
%   1. "DesignMultiSineExcitations" toolbox:
%
%           - Design of random phase multisine excitations with random harmonic 
%             grid (SISO systems only). 
%
%           - Design of full random orthogonal multisines for measuring MIMO systems. 
%
%           - Linear or logarithmic spacing of the excited harmonics.
%
%
%   2. "NL_Detection" toolbox:
%
%           - Detection, quantification, and classification of nonlinear distortions  
%             in FRF measurements using random phase multisine excitations 
%             (see "DesignMultiSineExcitations" toolbox). 
%
%           - Starts from the steady state response to (a) periodic input(s). 
%
%           - No noise transient (leakage) suppression. 
%
%           - Calculation of the sample means and sample covariances required by the 
%             sample maximum likelihood estimator of parametric transfer function models 
%             (see "MIMO_MaximumLikelihood" toolbox). 
%
%           - Fast (1 experiment with P >= 7 periods) and robust (M >= 7 experiments with 
%             P >= 2 periods) algorithms. 
%
%           - SISO systems only.
%
%
%   3. "LocalPolynomialMethod" toolbox:
%
%           - Detection and quantification of nonlinear distortions in FRF measurements using 
%             orthogonal multisine excitations (see "DesignMultiSineExcitations" toolbox). 
%             
%           - Can handle the transient response and multiple data sets.
%
%           - Nonparametric suppression of the plant and noise transients (leakage) errors 
%             in the FRF estimate (arbitrary and periodic signals).
%
%           - Calculation of the generalized sample means and sample covariances required 
%             by the sample maximum likelihood estimator of parametric transfer function 
%             models (see "MIMO_MaximumLikelihood" toolbox). 
%
%           - Fast (1 experiment with P >= 2 periods) and robust (M >= 2 experiments with 
%             P >= 2 periods) algorithms. 
%
%           - MIMO systems.
%
%
%   4. "MIMO_MaximumLikelihood" toolbox:
%
%           - Parametric plant transfer function modeling (common denominator model) starting
%             from the (generalized) sample means and sample (co-)variances of the noisy input 
%             noisy output DFT spectra (see "NL_Detection" and "LocalPolynomialMethod" toolboxes). 
%
%           - Errors-in-variables or generalized output error stochastic framework. 
%
%           - Calculation of the asymptotic covariance of the plant model parameters, the plant 
%             transfer function, the poles (+ resonance frequencies, damping ratios, and time 
%             constants), and the residue matrices (+ singular values, left and right singular vectors). 
%
%           - Discrete-time (z-domain), continuous-time (s-domain), and diffusion (sqrt(s)-domain)
%             transfer function models. 
%
%           - MIMO systems.
%
%
%   5. "MIMO_BoxJenkins" toolbox:
%
%           - Parametric plant and noise transfer function modeling (common denominator models) starting
%             from the known input and noisy output DFT spectra. The following model structures can be 
%             identified: ARMA, ARMAX, output error (OE), and Box-Jenkins (BJ). For systems operating in 
%             closed loop the feedback (controller) dynamics should be known; otherwise the estimates are biased. 
%
%           - Calculation of the asymptotic covariance of the plant and noise model parameters, the plant 
%             and noise transfer function, the corresponding poles (+ resonance frequencies, damping ratios, and time 
%             constants), and residue matrices (+ singular values, left and right singular vectors). 
%
%           - Discrete-time (z-domain), continuous-time (s-domain), and diffusion (sqrt(s)-domain)
%             plant and noise transfer function models. 
%
%           - MIMO systems.
%
%
%   Notes: 
%
%           1. For each toolbox example files are provided in the folder "examples". 
%
%           2. The auxiliary toolboxes "CovRootsResidues" and "MatrixArrayFunctions" contain functions used 
%              by the "MIMO_MaximumLikelihood" and  "MIMO_BoxJenkins" toolboxes.
%
%           3. For the SISO case of "MIMO_MaximumLikelihood" a commercial Matlab toolbox "fdident" with a graphical 
%              user interface is available. For more information see the website: 
%
%                               http://matlab.gamax.hu/english/products/fdident/ 
%   
%
%%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 21 October 2011 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%


