function [CovParam, Param] = CalcCovSingularVectors(CovParam, Param);
%
%   function [CovParam, Param] = CalcCovSingularVectors(CovParam, Param);
%
%
%   OUTPUT PARAMETERS
%
%       CovParam    =       struct{'res'} 
%                               CovParam.res    =   struct{'all', 'lsv', 'rsv'}
%                                                   CovParam.res.all    =   covariance matrix of (vec(Res))re; where ()re puts the real
%                                                                           and imaginary parts of the matrix on top of each other; vec()
%                                                                           stacks the columns of the matrix on top of each other; and Res
%                                                                           is a residue matrix of the ny x nu transfer function matrix G = B/A;
%                                                                           size: (2*ny*nu) x (2*ny*nu) x na
%                                                   CovParam.res.lsv    =   covariance left singular vectors of the residues [real(ur); imag(ur)];
%                                                                           size: 2*ny x 2*ny x min(ny, nu) x na
%                                                   CovParam.res.rsv    =   covariance right singular vectors of the residues [real(vr); imag(vr)];
%                                                                           size: 2*nu x 2*nu x min(ny, nu) x na
%
%       Param     =         struct{'res', 'poles', 'sv', 'lsv', 'rsv'}
%                               Param.res       =   residue matrices of the transfer function matrix G = B/A;
%                                                   size: ny x nu x na
%                                                   Note: the residues are sorted in the same order as ThePoles
%                               Param.poles     =   poles of the transfer function matrix G = B/A;
%                                                   size: na x 1
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
%       CovParam    =       struct{'res'} 
%                               CovParam.res    =   struct{'all'}
%                                                   CovParam.res.all    =   covariance matrix of (vec(Res))re; where ()re puts the real
%                                                                           and imaginary parts of the matrix on top of each other; vec()
%                                                                           stacks the columns of the matrix on top of each other; and Res
%                                                                           is a residue matrix of the ny x nu transfer function matrix G = B/A;
%                                                                           size: (2*ny*nu) x (2*ny*nu) x na
%
%       Param     =         struct{'res', 'poles', 'sv'}
%                               Param.res       =   residue matrices of the transfer function matrix G = B/A;
%                                                   size: ny x nu x na
%                                                   Note: the residues are sorted in the same order as ThePoles
%                               Param.poles     =   poles of the transfer function matrix G = B/A;
%                                                   size: na x 1
%                               Param.sv        =   singular values of the residue matrices;
%                                                   size: min(ny, nu) x na
%                               Param.lsv       =   left singular vectors of the residue matrices;
%                                                   size: ny x min(ny, nu) x na
%                               Param.rsv       =   right vectors of the residue matrices;
%                                                   size: nu x min(ny, nu) x na
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, May 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 23 March 2010
%

% initialisation variables
[ny, nu, na] = size(Param.res);
nsv = min([ny, nu]);            % number of singular values


% permutation matrix P from vec(A) to vec(A^H) where A is ny x nu
Enynu = zeros(ny, nu, ny, nu);
Enuny = zeros(nu, ny, nu, ny);
for ii=1:ny
	for jj=1:nu
		Enynu(ii, jj, ii, jj) = 1;
		Enuny(jj, ii, jj, ii) = 1;
	end % jj
end % ii

P = zeros(ny*nu, ny*nu);
for ii = 1:ny
	for jj=1:nu
		P = P + kron(squeeze(Enynu(:, :, ii, jj)), squeeze(Enuny(:, :, jj, ii)));
	end % jj
end % ii

% permutation matrix P1 from (vec(A))re to (vec(A^H))re where A is ny x nu
PRe = [real(P), -imag(P); imag(P), real(P)];
P1 = PRe * [eye(ny*nu), zeros(ny*nu); zeros(ny*nu), -eye(ny*nu)];

% initialisation covariance matrices
CovParam.res.lsv = zeros(2*ny, 2*ny, nsv, na);
CovParam.res.rsv = zeros(2*nu, 2*nu, nsv, na);

for ii = 1:na
    
    % residue matrix
     Aii = Param.res(:, :, ii);
     AiiRe = [real(Aii), -imag(Aii); imag(Aii), real(Aii)];

    for rr = 1:nsv
        
        ur = Param.lsv(:, rr, ii);                  % rr-th left singular vector
        vr = Param.rsv(:, rr, ii);                  % rr-th right singular vector
        [vrmax, irmax] = max(abs(vr));              % index first largest element vr
        phasor = exp(-sqrt(-1)*angle(vr(irmax)));
        ur = ur * phasor;
        vr = vr * phasor;
        Param.lsv(:, rr, ii) = ur;
        Param.rsv(:, rr, ii) = vr;
        sr = Param.sv(rr, ii);                      % rr-th singular value
        
        % Cr matrix depending on rr-th left and right singular vectors and rr-th singular value
		C11 = kron(vr.', eye(ny));
		C11 = [real(C11), -imag(C11); imag(C11), real(C11)];
		C12 = kron(conj(vr), ur);
		C12 = [real(ur); imag(ur)] * [real(C12); imag(C12)].';
		C21 = kron(ur.', eye(nu));
		C21 = [real(C21), -imag(C21); imag(C21), real(C21)];
		C22 = kron(conj(ur), vr);
		C22 = [real(vr); imag(vr)] * [real(C22); imag(C22)].';
		
		Cr = [C11 - C12; (C21 - C22) * P1]/sr;
        
        % Br-matrix
        Br = [eye(2*ny), -AiiRe/sr; -AiiRe.'/sr, eye(2*nu)];

        % Brplus is the pseudo-inverse of Br where imag(vr(irmax)) is removed
        Br(:, 2*ny+nu+irmax) = [];                   % impose imaginary part vr(irmax) to be zero
        [uBr, sBr, vBr] = svd(Br, 0);
        sBr = diag(sBr);
        Brplus = vBr * diag([1./sBr(1:end-1); 0]) * uBr.';
        
        % square root of CovParam.res.all
        [uu, ss, vv] = svd(squeeze(CovParam.res.all(:, :, ii)), 0);
        sqrtCov_res_all = uu * diag(diag(ss).^0.5);

        % covariance matrix of [ur; vr]
        BplusCr = Brplus * Cr * sqrtCov_res_all;
        Covurvr= BplusCr * BplusCr.';
%         Covurvr= BplusCr * CovParam.res.all(:, :, ii) * BplusCr.';
        CovParam.res.lsv(:, :, rr, ii) = Covurvr(1:2*ny, 1:2*ny);
        SelectIndex = [1:2*nu];
        SelectIndex(nu+irmax) = [];                   % remove imag(vr(irmax))
        CovParam.res.rsv(SelectIndex, SelectIndex, rr, ii) = Covurvr(2*ny+1:end, 2*ny+1:end);
       
    end % rr index singular value
    
end % ii index residue matrices
