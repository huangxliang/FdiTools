function [CovParam, Param] = OldCovResidues(Theta, CovTheta, Sel, ThePlane, Ts);
%
%   [CovParam, Param] = OldCovResidues(Theta, CovTheta, Sel, ThePlane, Ts);
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
%                                                   CovParam.res.all    =   covariance matrix of vec((Res)re); where ()re puts the real
%                                                                           and imaginary parts of the matrix on top of each other; vec()
%                                                                           stacks the columns of the matrix on top of each other; and Res
%                                                                           is a residue matrix of the ny x nu transfer function matrix G = B/A;
%                                                                           size: (2*ny*nu) x (2*ny*nu) x na
%                                                   CovParam.res.sv     =   variance singular values of the residues;
%                                                                           size: min(2*ny, nu) x na
%                                                   CovParam.re.lsv     =   cov((uk)re); where ()re puts the real and imaginary part of the kth
%                                                                           left singular vector uk on top of each other
%                                                                           size: (2*ny) x (2*ny) x nu x na
%                                                   CovParam.re.rsv     =   cov((vk)re); where ()re puts the real and imaginary part of the kth
%                                                                           right singular vector vk on top of each other
%                                                                           size: (2*nu) x (2*nu) x nu x na
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
%                                                   Note: the residues are sorted in the same order as ThePoles
%                               Param.poles     =   poles of the transfer function matrix G = B/A;
%                                                   size: na x 1
%                               Param.sv        =   singular values of the residue matrices;
%                                                   size: min(2*ny, nu) x na
%                               Param.lsv       =   left singular vectors of the residue matrices;
%                                                   size: ny x max(ny, nu) x na
%                               Param.rsv       =   left right vectors of the residue matrices;
%                                                   size: nu x max(ny, nu) x na
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
%   Rik Pintelon, February 2006
%   modified April 2006
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
[CovPoles, ThePoles] = CovRoots(Theta.A, CovTheta.A, Sel.A, ThePlane);
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
%   d vec((Res)re) / d (Poles)re;     d vec((Res)re) / d vec(B)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d Res / d Poles  and  d Res / d B w.r.t. ALL B-coefficients
% sorted in order of increasing magnitude of the poles
S_Res = SensResidues(Theta.B, Theta.A, ThePlane);

% d vec((Res)re) / d (Poles)re
S_Res_Poles = zeros(2*ny, nu, 2*na, na);
S_Res_Poles(:, :, 1:na, :) = [real(S_Res.poles); imag(S_Res.poles)];
S_Res_Poles(:, :, na+1:2*na, :) = [- imag(S_Res.poles); real(S_Res.poles)];
S_Res_Poles = reshape(S_Res_Poles, [2*ny*nu, 2*na, na]);

% d vec((Res)re) / d vec(B) w.r.t. ESTIMATED B-coefficients
% vec(B) is defined as follows: vecB = permute(B, [3, 1, 2]); vecB = vecB(:) 
SelB = permute(Sel.B, [3, 1, 2]);
SelB = SelB(:).';
S_Res.B = permute(S_Res.B, [1, 2, 5, 3, 4, 6]);                     % result: ny x nu x (nb+1) x ny x nu x na
S_Res_vecB = [real(S_Res.B); imag(S_Res.B)];
S_Res_vecB = reshape(S_Res_vecB, [2*ny*nu, (nb+1)*ny*nu, na]);      % sensitivity vec((Res)re) w.r.t. ALL vec(B)-coefficients
S_Res_vecB = S_Res_vecB(:, SelB == 1, :);                           % select the ESTIMATED vec(B)-coefficients


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix residues %
%       Cov((Res)re)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CovParam.res = struct('all', [], 'sv', [], 'lsv', [], 'rsv', []);
CovParam.res.all = zeros(2*ny*nu, 2*ny*nu, na);
for ii = 1:na
    CovParam.res.all(:, :, ii) = squeeze(S_Res_Poles(:, :, ii)) * CovPoles.all * squeeze(S_Res_Poles(:, :, ii)).' + ...
                                 squeeze(S_Res_vecB(:, :, ii)) * CovTheta.vecB * squeeze(S_Res_vecB(:, :, ii)).' + ...
                                 2 * herm(squeeze(S_Res_Poles(:, :, ii)) * CovPoles_vecB * squeeze(S_Res_vecB(:, :, ii)).');
end %ii


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation residues, their singular values, and the variance of the singular values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% in order of increasing magnitude of the poles
Param.res = CalcResidues(Theta.B, Theta.A, ThePlane);

% singular values and their variance
nsv = min(2*ny, nu);                % number of singular values
Param.sv = zeros(nsv, na);
CovParam.res.sv = zeros(nsv, na);
for ii = 1:na
    [uu, ss, vv] = svd([real(Param.res(:, :, ii)); imag(Param.res(:, :, ii))], 0); 
    Param.sv(:, ii) = diag(ss);
    for jj = 1:nsv
        TheMat = kron(vv(:, jj), uu(:, jj));
        CovParam.res.sv(jj, ii) = TheMat.' * squeeze(CovParam.res.all(:, :, ii)) * TheMat;
    end % jj
end % ii










