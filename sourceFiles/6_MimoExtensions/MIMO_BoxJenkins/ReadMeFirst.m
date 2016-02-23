% Read me first file: use of the routines for parametric plant and noise modeling of MIMO systems with nu inputs and ny outputs. 
% These functions use the "MatrixArrayFunctions" toolbox (MIMO_BoxJenkins), and the "CovRootsResidues" toolbox (MIMO_CR_bound).
% Examples illustrating the use of the functions can be found in the folder "examples". To run these examples one
% also needs the "LocalPolynomialMethod" and "MIMO_MaximumLikelihood" toolboxes.
%
%%   1. MAXIMUM LIKELIHOOD ESTIMATE PLANT AND NOISE MODEL PARAMETERS
%
%           [Theta, Sel, Cost, smax, smin, wscale] = MIMO_BoxJenkins(data, Sel, Theta0, ModelVar, IterVar) 
%
%       Starting from known input and noisy output DFT spectra, the function estimates the parameters of an ny x nu common denominator
%       plant model and an ny x ny common denominator noise model, and this for the following model structures:
%
%           1. ARMA
%           2. ARMAX
%           3. Output error (OE)
%           4. Box-Jenkins (BJ)
%
%       Discrete-time (z-domain), continuous-time (s-domain), as well as diffusion (sqrt(s)-domain) systems can be identified. 
%       If the plant operates in closed loop, then the feedback dynamics should be known; otherwise the estimates are biased. 
%
%       The structures Sel (definition free parameters), Theta0 (defintion parameters), ModelVar (definition model type), and 
%       IterVar (definition iteration parameters), as needed by the function MIMO_ML are easily generated using the function: 
%
%           [Sel, Theta0, ModelVar, IterVar] = MIMO_ML_DefaultValues(na, nb, nu, ny, PlantPlane, ModelStruct, Recip, Transient) 
%
%       The plant and noise transfer function matrices G(x, Theta) and H(s, Theta) corresponding to the model parameters Theta 
%       and the generalised frequency variable x (= j*omega for continuous-time; exp(-j*omega*Ts) for discrete-time, sqrt(j*omega) 
%       for diffusion) are calculated using the funtion:
%
%           PolyTrans = CalcPolyTrans(Theta, x);
%
%       Starting values for the plant and noise model parameters can be obtained with the following functions 
%
%           ThetaPlant = StartPlantModel(data, Sel, ModelVar, IterVar, FigNum)
%
%           Theta = StartNoiseModel(data, ThetaPlant, Sel, ModelVar, IterVar, FigNum)
%
%       References:
%
%           Pintelon R., and J. Schoukens (2006). Box-Jenkins identification revisited - Part I: theory, Automatica, 
%           vol. 42, no. 1, pp. 63-75.
%
%           Pintelon R., Y. Rolain, and J. Schoukens (2006). Box-Jenkins identification revisited - Part II: applications, 
%           Automatica, vol. 42, no. 1, pp. 77-84.
%
%           Pintelon, R., J. Schoukens, and P. Guillaume (2007). Box-Jenkins identification revisited - Part III: multivariable 
%           systems, Automatica, vol. 43, no. 5, pp. 868-875.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA). 
%
%
%%   2. CRAMER-RAO LOWER BOUND - ASYMPTOTIC COVARIANCE MATRIX
%
%           [CRbound, Theta, CovThetan, Thetan, Seln, wscale, TheCond] = MIMO_CR_bound(data, Sel, Theta, ModelVar) 
%
%       Calculates the Cramér-Rao lower bound of the physical model parameters (numerator and denominator coefficients, poles, 
%       residue matrices) of the estimated pnat and noise transfer function models. 
%
%       References:
%
%           Pintelon, R., P. Guillaume, and J. Schoukens (2007). Uncertainty calculation in (operational) modal analysis, 
%           Mechanical Systems and Signal Processing, vol. 21, no. 6, pp. 2359-2373.
%
%           Pintelon, R., J. Schoukens, and P. Guillaume (2007). Box-Jenkins identification revisited - Part III: multivariable 
%           systems, Automatica, vol. 43, no. 5, pp. 868-875.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA). 
%
%
%%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 5 November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% version 21 October 2011
%


