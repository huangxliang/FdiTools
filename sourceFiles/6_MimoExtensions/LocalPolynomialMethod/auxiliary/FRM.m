function [G, SUinv] = FRM(Z);
%
%       Calculates the frequency response matrix (FRM) starting from nu 
%       MIMO experiments using periodic signals. The inverse input matrix
%       is also calculated.
%
%
%   function [G, SUinv] = FRM(Z);
%
%
%	Output parameters
%
%
%       G       =   estimated frequency response matrix
%                   size ny x nu x F
%
%       SUinv	=   inverse of the input matrix inv(conj(U*U')) 
%                   size nu x nu x F
%
%
%	Input parameters
%
%       Z       =   matrix containing the output and input DFT spectra on top of each other of the nu MIMO experiments 
%                   size (ny+nu) x nu x F
%
%
% Rik Pintelon, October 21, 2009
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows, nu, F] = size(Z);
ny = rows - nu;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the frequency response matrix %
% and the inverse input matrix                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = zeros(ny, nu, F);
SUinv = zeros(nu, nu, F);
if nu == 1
    G(:, 1, :) = Z(1:ny, 1, :) ./ repmat(Z(end, 1, :), [ny, 1, 1]);
    SUinv(1, 1, :) = 1./squeeze(abs(Z(end, 1, :).^2));
else % nu > 1
    for kk = 1:F
        SUinv(:, :, kk) = inv(conj(squeeze(Z((ny+1):end, :, kk)) * squeeze(Z((ny+1):end, :, kk))'));
        G(:, :, kk) = Z(1:ny, :, kk) * inv(squeeze(Z((ny+1):end, :, kk)));
    end % kk
end % if one input
