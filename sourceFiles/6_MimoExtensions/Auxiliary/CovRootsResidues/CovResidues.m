function [CovParam, Param] = CovResidues(Theta, CovTheta, Sel, ThePlane, Ts);
%
%   [CovParam, Param] = CovResidues(Theta, CovTheta, Sel, ThePlane, Ts);
%
%       Calculates 
%                   1. the covariance matrix of the residue matrices of the transfer function matrix G = B/A
%                   2. the covariance matrix of the poles of G = B/A including the variance of the
%                      resonance frequencies, damping coefficients, and time constants
%                   3. the variance of the singular values of the residue matrices
%       where A is a polynomial of order na; and B is an ny x nu matrix polynomial of order nb
%       
%
%   OUTPUT PARAMETERS
%
%       CovParam    =       struct{'res', 'poles'} 
%                               CovParam.res    =   struct{'all', 'sv', 'lsv', 'rsv'}
%                                                   CovParam.res.all    =   covariance matrix of (vec(Res))re; where ()re puts the real
%                                                                           and imaginary parts of the matrix on top of each other; vec()
%                                                                           stacks the columns of the matrix on top of each other; and Res
%                                                                           is a residue matrix of the ny x nu transfer function matrix G = B/A;
%                                                                           size: (2*ny*nu) x (2*ny*nu) x na
%                                                   CovParam.res.sv     =   variance singular values of the residues;
%                                                                           size: min(ny, nu) x na
%                                                   CovParam.res.lsv    =   covariance left singular vectors of the residues [real(ur); imag(ur)];
%                                                                           size: 2*ny x 2*ny x min(ny, nu) x na
%                                                   CovParam.res.rsv    =   covariance right singular vectors of the residues [real(vr); imag(vr)];
%                                                                           size: 2*nu x 2*nu x min(ny, nu) x na
%                               CovParam.poles  =   struct{'root', 'all', 'damp', 'freq', 'time'}
%                                                   CovParam.poles.root =   cov((root)re); where ()re stacks the real and imaginary part of the root
%                                                                           on top of each other; ; size 2 x 2 x number of roots
%                                                   CovParam.poles.all  =   Croots.all = cov((roots)re); where the real and imaginary parts of
%                                                                           the vector roots are stacked on top of each other;
%                                                                           size 2*(number of roots) x 2*(number of roots)
%                                                   CovParam.poles.damp =   variance damping complex roots; entry is NaN for real roots
%                                                   CovParam.poles.freq =   variance frequency complex roots; entry is NaN for real roots
%                                                   CovParam.poles.time =   variance time constant real roots; entry is NaN for complex roots
%
%       Param     =         struct{'res', 'poles', 'sv'}
%                               Param.res       =   residue matrices of the transfer function matrix G = B/A;
%                                                   size: ny x nu x na
%                                                   Note: the residues are sorted in the same order as the poles 
%                               Param.poles     =   struct{'root', 'damp', 'freq', 'time}
%                                                   Param.poles.root    =   poles of the transfer function matrix G = B/A;
%                                                                           size na x 1
%                                                   Param.poles.damp    =   damping ratio of the poles;
%                                                                           size na x 1
%                                                   Param.poles.freq    =   resonance frequency of the poles;
%                                                                           size na x 1
%                                                   Param.poles.time    =   time constant of the real poles;
%                                                                           size na x 1
%                               Param.sv        =   singular values of the residue matrices;
%                                                   size: min(ny, nu) x na
%                               Param.lsv       =   left singular vectors of the residue matrices;
%                                                   size: ny x min(ny, nu) x na
%                               Param.rsv       =   right vectors of the residue matrices;
%                                                   size: nu x min(ny, nu) x na
%                               Note:               the largest element of the right singular vector is made positive real
%                                                   by multiplying the left and right singular vectors by exp(-j*fi) where
%                                                   fi is the phase of the largest element of the right singular vector
%
%
%   INPUT PARAMETERS
%
%       Theta           =   struct('A', 'B'} containing the numerator and denominator coefficients of the
%                           transfer function matrix G = B/A
%                               Theta.A =   denominator coefficients; size: 1 x (na+1)
%                               Theta.B =   numerator coefficients; size: ny x nu x (nb+1)
%
%       CovTheta        =   struct{'A', 'vecB', 'AvecB'} containing the covariance matrix of the estimated model parameters, 
%                           where vec(B) is calculated as: vecB = permute(B, [3, 1, 2]); vecB = vecB(:);
%                           CovTheta.A = cov(a.');                  size: nafree x nafree
%                           CovTheta.vecB = cov(vec(B));            size: nbfree x nbfree
%                           CovTheta.AvecB = covar(a.', vec(B));    size: nafree x nbfree
%                               nafree = number of estimated denominator coefficients
%                               nbfree = number of estimated numerator coefficients
%
%       Sel             =   struct('A', 'B') contains the position of the estimated coefficients
%                               Sel.A   =   size: 1 x (na+1)
%                                           Sel.A(ii) = 1 if Theta.A(ii) is estimated; otherwise zero
%                               Sel.B   =   size: ny x nu x (nb+1)
%                                           Sel.B(ii, jj, kk) = 1 if Theta.B(ii, jj, kk) is estimated; otherwise zero
%
%       ThePlane        =   domain of the polynomial
%					    	    's':	continuous-time;
%						        'w':	sqrt(s)-domain
%						        'z':	discrete-time;
%
%       Ts              =   sampling period in case of discrete-time models; optional; default value is one
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, February 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 5 October 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 4
    Ts = 1;
end % if
na = length(Sel.A) - 1;
[ny, nu, nb] = size(Sel.B);
nb = nb - 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix of the poles %
%       Cov((Poles)re)           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the poles are sorted in order of increasing magnitude
[CovPoles, ThePoles] = CovRoots(Theta.A, CovTheta.A, Sel.A, ThePlane, Ts);
CovParam.poles = CovPoles;
Param.poles = ThePoles;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Covariance between Poles and numerator coefficients   %
% Cov((Poles)re, vec(B)) = d(Poles)re/dA * Cov(A, vec(B) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sensitivity poles w.r.t. ALL denominator coefficients
% sorted in order of increasing magnitude of the poles
S_PolesA = SensRootsCoeff(Theta.A, ThePlane);

% sensitivity matrix of the poles w.r.t. the ESTIMATED A-coefficients
% d(Poles)re/dA
S_PolesA = S_PolesA(:, Sel.A == 1);
S_PolesA = [real(S_PolesA); imag(S_PolesA)];

% Cov((Poles)re, vec(B))
CovPoles_vecB = S_PolesA * CovTheta.AvecB;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity residues w.r.t. the poles and numerator coefficients %
%   d (vec(Res))re / d (Poles)re;     d (vec(Res))re / d vec(B)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d Res / d Poles  and  d Res / d B w.r.t. ALL B-coefficients
% sorted in order of increasing magnitude of the poles
% size S_Res.poles: ny x nu x na x na
% size S_Res.B: ny x nu x ny x nu x (nb+1) x na
S_Res = SensResidues(Theta.B, Theta.A, ThePlane);

% d vec(Res) / d Poles and d vec(Res) / d B
S_Res.poles = reshape(S_Res.poles, [ny*nu, na, na]);
S_Res.B = reshape(S_Res.B, [ny*nu, ny, nu, nb+1, na]);

% d (vec(Res))re / d (Poles)re
S_Res_Poles = zeros(2*ny*nu, 2*na, na);
S_Res_Poles(:, 1:na, :) = [real(S_Res.poles); imag(S_Res.poles)];
S_Res_Poles(:, na+1:2*na, :) = [- imag(S_Res.poles); real(S_Res.poles)];

% d (vec(Res))re / d vec(B) w.r.t. ESTIMATED B-coefficients
% vec(B) is defined as follows: vecB = permute(B, [3, 1, 2]); vecB = vecB(:) 
SelB = permute(Sel.B, [3, 1, 2]);
SelB = SelB(:).';
S_Res.B = permute(S_Res.B, [1, 4, 2, 3, 5]);                        % result: ny*nu x (nb+1) x ny x nu x na
S_Res_vecB = [real(S_Res.B); imag(S_Res.B)];                        % result: 2*ny*nu x (nb+1) x ny x nu x na
S_Res_vecB = reshape(S_Res_vecB, [2*ny*nu, (nb+1)*ny*nu, na]);      % sensitivity (vec(Res))re w.r.t. ALL vec(B)-coefficients
S_Res_vecB = S_Res_vecB(:, SelB == 1, :);                           % select the ESTIMATED vec(B)-coefficients


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix residues %
%     Cov((vec(Res))re)      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% square root covariance matrices
[uu, ss, vv] = svd([CovPoles.all, CovPoles_vecB; CovPoles_vecB.', CovTheta.vecB], 0);
sqrtCov = uu * diag(diag(ss).^0.5);

CovParam.res = struct('all', [], 'sv', []);
CovParam.res.all = zeros(2*ny*nu, 2*ny*nu, na);
for ii = 1:na
    J1 = [squeeze(S_Res_Poles(:, :, ii)), squeeze(S_Res_vecB(:, :, ii))] * sqrtCov;
    CovParam.res.all(:, :, ii) = J1 * J1.';
%     CovParam.res.all(:, :, ii) = squeeze(S_Res_Poles(:, :, ii)) * CovPoles.all * squeeze(S_Res_Poles(:, :, ii)).' + ...
%                                  squeeze(S_Res_vecB(:, :, ii)) * CovTheta.vecB * squeeze(S_Res_vecB(:, :, ii)).' + ...
%                                  2 * herm(squeeze(S_Res_Poles(:, :, ii)) * CovPoles_vecB * squeeze(S_Res_vecB(:, :, ii)).');
end %ii


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation residues, their singular vectors, singular values, and the variance of the singular values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% in order of increasing magnitude of the poles
Param.res = CalcResidues(Theta.B, Theta.A, ThePlane);

% singular values and their variance
nsv = min(ny, nu);                % number of singular values
Param.sv = zeros(nsv, na);
Param.lsv = zeros(ny, nsv, na);
Param.rsv = zeros(nu, nsv, na);
CovParam.res.sv = zeros(nsv, na);
for ii = 1:na
    if ny < nu
       [vv, ss, uu] = svd(Param.res(:, :, ii)', 0); 
    else
       [uu, ss, vv] = svd(Param.res(:, :, ii), 0); 
    end % if
    Param.sv(:, ii) = diag(ss);
    Param.lsv(:, :, ii) = uu;
    Param.rsv(:, :, ii) = vv;
    [u1, s1, v1] = svd(squeeze(CovParam.res.all(:, :, ii)), 0);
    sqrtCov_res_all = u1 * diag(diag(s1).^0.5);
    for jj = 1:nsv
        TheMat = kron(conj(vv(:, jj)), uu(:, jj));
        TheMat = sqrtCov_res_all.' * [real(TheMat); imag(TheMat)];
        CovParam.res.sv(jj, ii) = TheMat.' * TheMat;
%         CovParam.res.sv(jj, ii) = TheMat.' * squeeze(CovParam.res.all(:, :, ii)) * TheMat;
    end % jj
end % ii index residue matrices


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation covariance matrix left and right singular vectors of the residue matrices %
% Note: the largest element of the right singular vector is made positive real by       %
%       by multiplying the left and right singular vectors by exp(-j*fi), with fi the   %
%       phase of the largest element                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[CovParam, Param] = CalcCovSingularVectors(CovParam, Param);

