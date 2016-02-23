%% 1. 
% SET DEFAULT VALUES:
%
%Sel =	selection of the estimated and the known coefficients
%Sel = struct('A', [],'B', [], 'Ig', [])
%Theta0	= starting value plant, noise, and initial conditions parameters
%Theta0 = struct('A', [], 'B', [], 'Ig', [])
%ModelVar =	contains the information about the model to be identified
%ModelVar = struct(Transient=1, PlantPlane=s/z, Struct=EIV/OE, Reciprocal=0)
%IterVar = contains the information about the minimization procedure
%IterVar = struct(LM=1,MaxIter=100,TolParam=1e-6,TolCost=1e-10,TraceOn=1,NormJacob=1)

%% 2.
% CALCULATE NON-PARA FRM
%
%CZ = sample covariance matrices of the sample mean of the Fourier coefficients Z
%   CZ.n = noise covariance matrix of Z as complex numbers
%   CZ.NL = total (noise + NL distortion) covariance matrix of Z
%   CZ.S = covariance matrix of the stochastic nonlinear distortions in Z

%Z = sample mean Fourier coefficients
%freq = frequency of the excited harmonics; size 1 x F
%G = estimated frequency response matrix; size ny x nu x F
%CvecG = covariance matrices of vec(G); size (ny*nu) x (ny*nu) x F

%   CvecG.n = noise covariance matrix of vec(G)
%   CvecG.NL = total (noise + NL distortion) covariance matrix of vec(G)
%   CvecG.NLcheck = total (noise + nonlinear distortion) covariance matrix of vec(G) calculated
%                   via the classical averaging over the realisations 
%                   (no averaging over neighbouring excited frequencies)
%   CvecG.S = covariance matrix of the stochastic NL distortions in vec(G)

%dof = actual degrees of freedom of the input-output and FRF covariance estimates 
%    = number of equivalent independent experiments - 1
%CL = ± CL is the correlation length (over the frequency) of 
%       the sample mean and/or sample covariance
% => G is FRF: mag*exp(phs*j*pi/180)

%% 3.
% CALCULATE PARA ML
%
%Theta = estimated value plant, noise, and initial conditions parameters
%Cost  = value of the cost function in the last iteration step
%smax  = largest singular value of the Jacobian matrix
%smin  = smallest singular value of the Jacobian matrix
%wscale	= angular frequency scaling
% => GML is H_ml

%% 4.
% CRAMER-RAO LOWER BOUND
%
%CRbound = Cramer-Rao bound of the estimated plant model parameters
%CRbound = struct('A', [], 'AvecB', [], 'vecB', [], 'vecG', [], 'res', [], 'poles')


%% 5.
% WHITENESS TEST RESIDUALS