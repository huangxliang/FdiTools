function W = MIMO_Calc_IQML_Weight(data, PolyTrans);
%
%       Calculates the noise weighting of the linear least squares cost function 
%
%   function W = MIMO_Calc_IQML_Weight(data, PolyTrans);
%
%
%   Output
%
%       W           =   weighting matrix for the linear least squares cost function
%                           1 MIMO experiment:   ny x ny x F 
%                           nu MIMO experiments: ny x ny x nu x F 
%
%   Input
%
%		data		=	structure containing the non-parametric data required for the identification
%							data.CY         =	(sample) noise covariance matrix of Y 
%                                               	1 MIMO experiment:                  ny x ny x F 
%                                                	nexp MIMO experiments (nexp >= 1):	ny x ny x nexp x F 
%                           data.CU         =   (sample) noise covariance matrix of U  
%                                                   1 MIMO experiment:   nu x nu x F 
%                                                	nexp MIMO experiments (nexp >= 1):  nu x nu x nexp x F 
%                           data.CYU        =   (sample) noise covariance matrix of U 
%                                                	1 MIMO experiment:   ny x nu x F 
%                                               	nexp MIMO experiments (nexp >= 1):  ny x nu x nexp x F 
%
%                           data.G          =   frequency response matrix, size ny x nu x F 
%
%		PolyTrans	=	structure containing the polynomials and transfer functions evaluated in x
%							PolyTrans.A		=	denominator polynomial plant transfer function evaluated in x.Plant, size 1 x F 
%							PolyTrans.G		=	plant transfer matrix evaluated in x.Plant, size ny x nu x F
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 26 November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 24 October 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the weighting matrix W %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine the size of the covariance matrices
NumberDim = length(size(data.CY));

[ny, nu, F] = size(PolyTrans.G);
if NumberDim == 4       % nexp MIMO experiments
    
    nexp = size(data.CY, 3);
    W = zeros(ny, ny, nexp, F);
    dummy0 = zeros(nu, nu, F);
    dummy1 = zeros(ny, nu, F);
    dummy2 = zeros(ny, ny, 1, F);
    for ee = 1:nexp
        dummy0(:,:,:) = data.CU(:,:,ee,:);
        dummy1(:,:,:) = data.CYU(:,:,ee,:);
        dummy2(:,:,1,:) = Mat_Mult(PolyTrans.G, Mat_Mult(dummy0, Conj_Trans(PolyTrans.G))) ...
                          - 2*herm(Mat_Mult(dummy1, Conj_Trans(PolyTrans.G)));
        W(:,:,ee,:) = data.CY(:,:,ee,:) + dummy2;           
    end % ee, nexp MIMO experiments
    W = W .* repmat(reshape(abs(PolyTrans.A.^2), [1,1,1,F]), [ny,ny,nexp,1]);
    
elseif NumberDim == 3	% one MIMO experiment
    
    W = data.CY + Mat_Mult(PolyTrans.G, Mat_Mult(data.CU, Conj_Trans(PolyTrans.G))) ...
        - 2*herm(Mat_Mult(data.CYU, Conj_Trans(PolyTrans.G))); 
    W = W .* repmat(reshape(abs(PolyTrans.A.^2), [1,1,F]), [ny,ny,1]);
    
end % if
