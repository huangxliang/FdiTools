% Read me first file: use of the routines for parametric plant modeling of MIMO systems with nu inputs and ny outputs. 
% These functions use the "MatrixArrayFunctions" toolbox (MIMO_ML, MIMO_WTLS, MIMO_WGTLS, MIMO_IQML, MIMO_BTLS, MIMO_ML_CRbound), 
% and the "CovRootsResidues" toolbox (MIMO_ML_CRbound).
% Examples illustrating the use of the functions can be found in the folder "examples". To run these examples one
% also needs the "LocalPolynomialMethod" and "DesignMultiSineExcitations" toolboxes.
%
%%   1. (SAMPLE) MAXIMUM LIKELIHOOD ESTIMATE PLANT MODEL PARAMETERS
%
%           [Theta, Cost, smax, smin, wscale] = MIMO_ML(data, Sel, Theta0, ModelVar, IterVar);
%
%       Starting from noisy input and noisy output DFT spectra (= the sample means) and their sample covariance based on dof
%       degrees of freedom, the function estimates the model parameters Theta of an ny x nu common denominator plant transfer 
%       function of order nb over na. Discrete-time (z-domain), continuous-time (s-domain), as well as diffusion (sqrt(s)-domain)
%       systems can be identified. The sample means and sample covariances are the output of the following functions:
%
%           1. ArbLocalPolyAnal:    1 MIMO experiment with random inputs (can be the concatenation of several data records)
%           2. FastLocalPolyAnal:	1 MIMO experiment with periodic inputs 
%                                   (random phase multisines, P > 1 periods) 
%           3. RobustLocalPolyAnal:	nu MIMO experiments with periodic inputs 
%                                   (random phase multisines, M >=1 realisations, P > 1 periods) 
%
%       Note: the MIMO_ML function can handle nexp >= 1 experiments with periodic signals.
%       In the absence of modeling errors the expected value and the variance of the cost function are given by 
%
%           E{Cost} = dof/(dof-ny)*ny*(F-ntheta/2)	  and       var_min	<= var{Cost} <= 3*var_min 
%
%       where
%
%           var_min = dof^3/((dof-ny)^2*(dof-ny-1))*ny*(F-ntheta/2)
%
%       The upper bound of the variance is reached for 1 MIMO experiment with random or periodic excitations. 
%       For nu MIMO experiments with periodic excitations (or 1 MIMO experiment with periodic excitation and averaging
%       over the periods only) the variance is in between the two bounds. 
%
%       The structures Sel (definition free parameters), Theta0 (defintion parameters), ModelVar (definition model type), and 
%       IterVar (definition iteration parameters), as needed by the function MIMO_ML are easily generated using the function: 
%
%           [Sel, Theta0, ModelVar, IterVar] = MIMO_ML_DefaultValues(na, nb, nu, ny, PlantPlane, ModelStruct, Recip, Transient) 
%
%       The transfer function matrix G(x, Theta) corresponding to the model parameters Theta and the generalised frequency 
%       variable x (= j*omega for continuous-time; exp(-j*omega*Ts) for discrete-time, sqrt(j*omega) for diffusion) is calculated
%       using the funtion:
%
%           PolyTrans = MIMO_ML_CalcPolyTrans(Theta, x);
%
%       To perform the whitness test on the residuals one should take into account the correlation ± CL over the frequency 
%       of the sample mean and sample covariance: the autocorrelation should be performed on residuals in steps of CL+1. 
%       over the EXCITED frequencies. This is automatically done in the function: 
%
%           [Auto_Corr, Lags, Conf_Bound, Fraction] = WhitenessTestResiduals(G, CvecG, Gest, CL, dof, PlotFig); 
%
%       References:
%
%           Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and FRF models 
%           for multivariable systems - Part II: extensions, applications, Mechanical Systems and Signal Processing, vol. 24, 
%           no. 3, pp. 596-616.
%
%           Pintelon, R., G. Vandersteen, J. Schoukens, and Y. Rolain (2011). Improved (non-)parametric identification of dynamic 
%           systems excited by periodic signals - The multivariate case, Mechanical Systems and Signal Processing, vol. 25, no. 8, 
%           pp. 2892-2922.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA). 
%
%
%%   2. CRAMER-RAO LOWER BOUND - ASYMPTOTIC COVARIANCE MATRIX
%
%           [CRbound, G, Theta, CovThetan, Thetan, Seln, wscale, TheCond] = MIMO_ML_CRbound(data, Sel, Theta, ModelVar); 
%
%       Calculates the Cramér-Rao lower bound of the physical model parameters (numerator and denominator coefficients, poles, 
%       residue matrices) and the estimated transfer function model for the following cases: 
%
%           1. 1 MIMO experiment with random excitations (the input or the power spectrum of the input is provided); 
%              or the concatenation of several data records (see the ArbLocalPolyAnal function).  
%           2. 1 MIMO experiment with periodic excitations
%           3. nexp >=1 MIMO experiments with periodic excitations
%
%       The covariance matrix of the SML estimate Theta_SML(Zm) using the sample mean Zm and sample covariance CZm (result of
%       the functions ArbLocalPolyAnal, or FastLocalPolyAnal) is related to the covariance matrix of the ML estimate Theta_ML(Z) 
%       using the original data set Z and the known covariance CZ0 as 
%
%           Cov(Theta_SML(Zm)) = dof*(dof-ny)/((dof-ny+1)*(dof-ny-1)) * Cov(Theta_ML(Z)) 
%
%       where Cov(Theta_ML(Z)) = CRbound.Theta(Z, CZ0) (= CRbound.Theta with (Z, CZ0) as input data). For the result of 
%       RobustLocalPolyAnal (or FastLocalPolyAnal with averaging over the periods only) Cov(Theta_ML(Z)) in the formula above is 
%       replaced by CRbound.Theta(Zm, CZm0) with CZm0 the true covariance of Zm. Since in practise the true covariance CZ0 (or 
%       CZm0) is unknown, Cov(Theta_ML(Z)) can only be calculated in an estimate CZ (sample cov.) of the true covariance CZ0. 
%       The formula above is then replaced by 
%
%           Cov(Theta_SML(Zm)) = dof^2/((dof-ny+1)*(dof-ny-1)) * CRbound(Z, CZ) 
%
%       where CRbound(Z, CZ) is replaced by CRbound(Zm, CZm) for RobustLocalPolyAnal  (or FastLocalPolyAnal with averaging 
%       over the periods only) 
%
%       References:
%
%           Pintelon, R., P. Guillaume, and J. Schoukens (2007). Uncertainty calculation in (operational) modal analysis, 
%           Mechanical Systems and Signal Processing, vol. 21, no. 6, pp. 2359-2373.
%
%           Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and FRF models 
%           for multivariable systems - Part II: extensions, applications, Mechanical Systems and Signal Processing, vol. 24, 
%           no. 3, pp. 596-616.
%
%           Pintelon, R., G. Vandersteen, J. Schoukens, and Y. Rolain (2011). Improved (non-)parametric identification of dynamic 
%           systems excited by periodic signals - The multivariate case, Mechanical Systems and Signal Processing, vol. 25, no. 8, 
%           pp. 2892-2922.
%
%           Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%           IEEE Press-Wiley, Piscataway (USA). 
%
%
%%  STARTING VALUES: SAMPLE WEIGHTED (GENERALIZED) TOTAL LEAST SQUARES, BOOTSTRAPPED TOTAL LEAST SQUARES 
%
%           [Theta, smax, smin, wscale] = MIMO_WTLS(data, Sel, ModelVar)
%
%       Weighted total least squares (WTLS) estimate of the common denominator plant model parameters. The routine does not
%       require any noise covariance information. In general the estimates are inconsistent. 
%
%           [Theta, smax, smin, wscale] = MIMO_WGTLS(data, Sel, ModelVar)
%
%       Weighted generalized total least squares (WGTLS) esimate of the common denominator plant model parameters. The routine requires
%       the (sample) noise covariance information. The estimates are consistent.
%
%           [Theta, Cost, smax, smin, wscale] = MIMO_IQML(data, Sel, Theta0, ModelVar, IterVar) 
%
%       Iterative quadratic maximum likelihood (IQML) estimate of the common denominator plant model parameters. The routine  
%       requires starting values and the (sample) noise covariance information. In general the estimates are inconsistent. 
%
%           [Theta, Cost, smax, smin, wscale] = MIMO_BTLS(data, Sel, Theta0, ModelVar, IterVar) 
%
%       Bootstrapped total least squares (BTLS) estimate of the common denominator plant model parameters. The routine requires 
%       starting values and the (sample) noise covariance information. The estimates are consistent. 
%
%       References:
%
%           Pintelon R., P. Guillaume, G. Vandersteen and Y. Rolain (1998). Analyses, development and applications of TLS 
%           algorithms in frequency-Domain System Identification, SIAM J. Matrix Anal. Appl., vol. 19, no. 4, pp. 983-1004.
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
% version 24 October 2011
%


