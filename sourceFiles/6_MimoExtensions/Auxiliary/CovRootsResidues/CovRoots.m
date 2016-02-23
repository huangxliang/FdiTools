function [CovTheRoots, TheRoots] = CovRoots(Coeff, CovCoeff, Sel, ThePlane, Ts);
%
%   [CovTheRoots, TheRoots] = CovRoots(Coeff, CovCoeff, Sel, ThePlane, Ts);
%
%
%       Calculates the covariance matrix of the roots of a polynomial defined by Coeff
%
%               Coeff(1) + Coeff(2)*x + Coeff(3) * x^2 + ... + Coeff(nx+1) x^nx
%
%       s-, sqrt(s)-domain  :   x = s, sqrt(s)
%       z-domain            :   x = z^-1
%
%
%   OUTPUT PARAMETER
%
%       CovTheRoots     =  structure containing the covariance of the roots struct{'root', 'all', 'damp', 'freq', 'time'}
%                               CovTheRoots.root    =   cov((root)re); where re stacks the real and imaginary part of the root
%                                                       on top of each other; ; size 2 x 2 x number of roots
%                               CovTheRoots.all     =   Croots.all = cov((roots)re); where the real and imaginary parts of the vector roots are
%                                                       stacked on top of each other; size 2*(number of roots) x 2*(number of roots)
%                               CovTheRoots.damp    =   variance of the damping of the complex roots; entry is NaN for real roots
%                               CovTheRoots.freq    =   variance of the frequency of the complex roots; entry is NaN for real roots
%                               CovTheRoots.time    =   variance of the time constant of the real roots; entry is NaN for complex roots
%                           Notes:
%                                   - in s, sqrt(s)-domain CovTheRoots is the covariance matrix of the frequency normalised roots
%                                   - everything is calculated in order of the vector TheRoots 
%
%       TheRoots        =   structure containing the roots of the polynomial, and related quantities such as damping, frequency,
%                           time constant struct('root', 'damp', 'freq', 'time')
%                               TheRoots.root       =   roots of the polynomial
%                               TheRoots.damp       =   damping of the complex roots; entry is NaN for real roots
%                               TheRoots.freq       =   frequency of the complex roots; entry is NaN for real roots
%                               TheRoots.time       =   time constant of the real roots; entry is NaN for complex roots
%                           Note: the roots are sorted in increasing order of magnitude
%
%
%   INPUT PARAMETERS
%
%       Coeff           =   coefficients of a polynomial in increasing powers of s, sqrt(s), or z^-1; size 1 x (1 + order polynomial)
%
%       CovCoeff        =   covariance matrix of the estimated polynomial coefficients; size (number of estimated coeff.) x (number of estimated coeff.)
%
%       Sel             =   Sel(i) = 1 if Coeff(i) is estimated
%                           Sel(i) = 0 if Coeff(i) is known
%
%       ThePlane        =   domain of the polynomial
%					    	    's':	continuous-time;
%						        'w':	sqrt(s)-domain
%						        'z':	discrete-time;
%
%       Ts              =   sampling period for polynomials in z^-1; otherwise may be empty
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, January 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 6 October 2011
%


%%%%%%%%%%%%%%%%%%
% initialisation %
%%%%%%%%%%%%%%%%%%

% order polynomial
nn = length(Coeff) - 1;

% number of estimated coefficients
ntheta = sum(Sel);

% initialise sample period if not given
if nargin == 4
    Ts = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity of the roots w.r.t. the estimated polynomial coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sensitivity matrix of the roots w.r.t. all polynomial coefficients
[S_RootsCoeff, RootsPoly] = SensRootsCoeff(Coeff, ThePlane);
TheRoots.root = RootsPoly;

% sensitivity matrix of the roots w.r.t. the estimated polynomial coefficients
S_RootsCoeff = S_RootsCoeff(:, Sel == 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the covariance matrices of the roots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% square root of the covariance matrix of the coefficients
[uu, ss, vv] = svd(CovCoeff, 0);
sqrtCovCoeff = uu * diag(diag(ss).^0.5);

% Covariance matrix of all the roots
CovTheRoots.all = zeros(2*nn, ntheta);
Sall = [real(S_RootsCoeff); imag(S_RootsCoeff)] * sqrtCovCoeff;
CovTheRoots.all = Sall * Sall.';

% Covariance matrix of each root seperately
CovTheRoots.root = zeros(2, 2, nn);
for ii = 1:nn
    Sroot = [real(S_RootsCoeff(ii, :)); imag(S_RootsCoeff(ii, :))] * sqrtCovCoeff;
    CovTheRoots.root(:, :, ii) = Sroot * Sroot.';    
end % ii


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Related quantities:                                 %
%   - damping                                         %
%   - frequency                                       %
%   - time constant                                   %
% In sqrt(s)-domain the roots are first transformed   %
% transformed using s = w^2.                          %
% In z-domain the roots are first transformed using   %
% the impulse invariant transformation s = log(z)/Ts. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ThePlane
    
    case 's'
        NewRoots = TheRoots.root;
        CovNewRoots = CovTheRoots.root;
        
    case 'w'
        % s = lambda^2 transformation
        NewRoots = TheRoots.root.^2;
        CovNewRoots = zeros(2, 2, nn);
        % transformation covariance matrix of each root
        for ii = 1:nn
            [uu, ss, vv] = svd(squeeze(CovTheRoots.root(:,:,ii)), 0);
            sqrtCovTheRoots = uu * diag(diag(ss).^0.5);
            Sp_z = 2*[real(TheRoots.root(ii)), -imag(TheRoots.root(ii)); imag(TheRoots.root(ii)), ...
                    real(TheRoots.root(ii))] * sqrtCovTheRoots;
            CovNewRoots(:,:,ii) = Sp_z * Sp_z.';
        end % ii
        
    case 'z'
        % impulse invariant transformation s = 1/Ts * log(z)
        NewRoots = log(TheRoots.root)/Ts;
        CovNewRoots = zeros(2, 2, nn);
        % transformation covariance matrix of each root
        for ii = 1:nn
            [uu, ss, vv] = svd(squeeze(CovTheRoots.root(:,:,ii)), 0);
            sqrtCovTheRoots = uu * diag(diag(ss).^0.5);
            Sp_z = [real(TheRoots.root(ii)), imag(TheRoots.root(ii)); -imag(TheRoots.root(ii)), ...
                    real(TheRoots.root(ii))] * sqrtCovTheRoots /(Ts*abs(TheRoots.root(ii)^2));
            CovNewRoots(:,:,ii) = Sp_z * Sp_z.';
        end % ii
     
end % switch

% find real and complex roots and initialize the variables
IndexReal = find(imag(NewRoots) == 0);                          % real roots
IndexCompl = find(imag(NewRoots) ~= 0);                         % complex roots
TheRoots.damp = zeros(size(TheRoots.root)) + NaN;
TheRoots.freq = zeros(size(TheRoots.root)) + NaN;
TheRoots.time = zeros(size(TheRoots.root)) + NaN;

% calculate damping, frequency, and time constant
TheRoots.freq(IndexCompl) = abs(NewRoots(IndexCompl))/(2*pi);
TheRoots.damp(IndexCompl) = - real(NewRoots(IndexCompl))./(2*pi*TheRoots.freq(IndexCompl));
TheRoots.time(IndexReal) = - 1./NewRoots(IndexReal);

% variance damping
CovTheRoots.damp = zeros(size(TheRoots.root)) + NaN;
ncompl = length(IndexCompl);
for ii = 1:ncompl
    [uu, ss, vv] = svd(squeeze(CovNewRoots(:,:,IndexCompl(ii))), 0);
    sqrtCovNewRoots = uu * diag(diag(ss).^0.5);
    VV = [- imag(NewRoots(IndexCompl(ii))), real(NewRoots(IndexCompl(ii)))] * sqrtCovNewRoots;
    CovTheRoots.damp(IndexCompl(ii)) = (VV * VV.') * ...
                                       (imag(NewRoots(IndexCompl(ii)))^2 / (2*pi*TheRoots.freq(IndexCompl(ii)))^6);
end
 
% variance frequency
CovTheRoots.freq = zeros(size(TheRoots.root)) + NaN;
ncompl = length(IndexCompl);
for ii = 1:ncompl
    [uu, ss, vv] = svd(squeeze(CovNewRoots(:,:,IndexCompl(ii))), 0);
    sqrtCovNewRoots = uu * diag(diag(ss).^0.5);
    VV = [real(NewRoots(IndexCompl(ii))), imag(NewRoots(IndexCompl(ii)))] * sqrtCovNewRoots / (2*pi)^2;
    CovTheRoots.freq(IndexCompl(ii)) = (VV * VV.') / TheRoots.freq(IndexCompl(ii))^2;
end

% variance time constant
CovTheRoots.time = zeros(size(TheRoots.root)) + NaN;
nreal = length(IndexReal);
for ii = 1:nreal
    CovTheRoots.time(IndexReal(ii)) = CovNewRoots(1,1,IndexReal(ii)) / NewRoots(IndexReal(ii))^4 ;
end
