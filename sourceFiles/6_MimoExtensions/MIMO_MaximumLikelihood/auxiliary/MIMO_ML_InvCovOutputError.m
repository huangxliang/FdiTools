function PolyTrans = MIMO_ML_InvCovOutputError(data, PolyTrans);
%
%           Calculates a hermitian symmetric square root of the inverse of 
%           the covariance matrix of the output error Y-G*U.
%
%   function PolyTrans = MIMO_ML_InvCovOutputError(data, PolyTrans);
%
%
%   Output
%
%		PolyTrans	=	structure containing the polynomials and transfer functions evaluated in x
%							PolyTrans.A             =	denominator polynomial plant transfer function evaluated in x.Plant 
%                                                       size 1 x F 
%							PolyTrans.G             =	plant transfer matrix evaluated in x.Plant
%                                                       size ny x nu x F 
%							PolyTrans.Tg            =	plant transient term evaluated in x.Plant
%                                                       size ny x F 
%                           PolyTrans.sqrtCEinv     =   hermitian symmetric square root of the inverse of the covariance of the 
%                                                       output error (Cov(NY-G*NU)) 
%
%   Input
%
%		data		=	structure containing the non-parametric covariance data
%							data.CY                 =	(sample) noise covariance matrix of Y, size: ny x ny x F 
%                           data.CU                 =   (sample) noise covariance matrix of U, size: nu x nu x F 
%                           data.CYU                =   (sample) noise covariance matrix of U, size: ny x nu x F 
%
%		PolyTrans	=	structure containing the polynomials and transfer functions evaluated in x
%							PolyTrans.A             =	denominator polynomial plant transfer function evaluated in x.Plant 
%                                                       size 1 x F 
%							PolyTrans.G             =	plant transfer matrix evaluated in x.Plant
%                                                       size ny x nu x F 
%							PolyTrans.Tg            =	plant transient term evaluated in x.Plant
%                                                       size ny x F 
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 1 December 2009
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = size(data.CY, 3);                       % number of frequencies
ny = size(data.CY, 1);                      % number of outputs

PolyTrans.sqrtCEinv = zeros(ny, ny, F);     % hermitian square root of CEinv 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast calculation of the covariance Cov(NY-G*NU) of the equation error (NY-G*NU). The lines     %
% below are equivalent with                                                                      % 
% for kk = 1:F                                                                                   % 
%	CE = data.CY(:,:,kk) + PolyTrans.G(:,:,kk) * data.CU(:,:,kk) * PolyTrans.G(:,:,kk)' ...      % 
%         - 2*herm(data.CYU(:,:,kk) * PolyTrans.G(:,:,kk)');                                     % 
% end                                                                                            % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PolyTrans.sqrtCEinv is used as intermediate variable for CE = Cov(NY-G*NU) 
PolyTrans.sqrtCEinv = data.CY + Mat_Mult(Mat_Mult(PolyTrans.G, data.CU), Conj_Trans(PolyTrans.G)) ...
                      - 2*herm(Mat_Mult(data.CYU, Conj_Trans(PolyTrans.G)));

                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of a hermitian symmetric square root of CEinv = Cov(Ny-G*NU) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ny == 1
    PolyTrans.sqrtCEinv = PolyTrans.sqrtCEinv.^(-0.5);
else % more than one output
    % relative upper limit under which the singular values are set to zero
    % in the pseudo-inverse of the square root of the covariance matrix
    delta = 1e-6; 
    % pseudo-inverse covariance matrix output residuals
    PolyTrans.sqrtCEinv = Sqrt_Inv(PolyTrans.sqrtCEinv, delta);
end % if one output

