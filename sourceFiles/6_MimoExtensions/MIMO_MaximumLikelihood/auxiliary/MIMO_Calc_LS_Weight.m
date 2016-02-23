function W = MIMO_Calc_LS_Weight(data, G);
%
%       Calculates the noise weighting of the linear least squares cost function 
%
%   function W = MIMO_Calc_LS_Weight(data, G);
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
%                                               	1 MIMO experiment:   ny x ny x F 
%                                                	nu MIMO experiments: ny x ny x nu x F 
%                           data.CU         =   (sample) noise covariance matrix of U  
%                                                   1 MIMO experiment:   nu x nu x F 
%                                                	nu MIMO experiments: nu x nu x nu x F 
%                           data.CYU        =   (sample) noise covariance matrix of U 
%                                                	1 MIMO experiment:   ny x nu x F 
%                                               	nu MIMO experiments: ny x nu x nu x F 
%
%                           data.G          =   frequency response matrix, size ny x nu x F 
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 26 November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine the number of MIMO experiments %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumberDim = length(size(data.CY));          % number of matrix dimensions
if NumberDim == 3
    NumberExp = 1;                          % number of MIMO experiments
elseif NumberDim == 4
    NumberExp = size(data.CY, 3);           % number of MIMO experiments
end % if


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the weighting matrix W %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NumberExp > 1
    
    [ny, nu, F] = size(G);
    W = zeros(ny, ny, nu, F);
    dummy1 = zeros(ny, nu, F);
    dummy2 = zeros(ny, ny, 1, F);
    for ee = 1:NumberExp
        dummy1(:,:,:) = data.CYU(:,:,ee,:);
        dummy2(:,:,1,:) = Mat_Mult(G, Mat_Mult(squeeze(data.CU(:,:,ee,:)), Conj_Trans(G))) ...
                          - 2*herm(Mat_Mult(dummy1, Conj_Trans(G)));
        W(:,:,ee,:) = data.CY(:,:,ee,:) + dummy2;           
    end % ee, MIMO experiments
    
else % one MIMO experiment
    
    W = data.CY + Mat_Mult(G, Mat_Mult(data.CU, Conj_Trans(G))) - 2*herm(Mat_Mult(data.CYU, Conj_Trans(G)));       
    
end % if more than one MIMO experiment
