function CvecG = Cov_FRM(G, SUInv, CZ);
%
%       Calculates the covariance matrix of vec(G) starting from
%       the covariance matrix of the input-output spectra
%
%   function CvecG = Cov_FRM(G, SUInv, CZ);
%
%
%   Output parameters
%
%       CvecG   =   covariance matrix of vec(G), taking into account that the nu experiments
%                   are independent with equal covariance matrix
%                   size ny*nu x ny*nu x F 
%
%
%   Input parameters
%
%       G       =   frequency response matrix
%                   size ny x nu x F
%
%       SUInv   =   inverse of the input matrix inv(conj(U * U'))
%                   size nu x nu x F
%
%       CZ      =   covariance matrix of the output and input spectra on top of each other of the nu experiments 
%                   size (ny+nu) x (ny+nu) x nu x F
%
% Rik Pintelon, October 22, 2009
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny, nu, F] = size(G);
CvecG = zeros(ny*nu, ny*nu, F);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of cov(vec(G)) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% average the covariance matrix over the nu experiments
CZ = squeeze(mean(CZ, 3));

if nu*ny == 1       % SISO case
    
    CvecG(1, 1, :) = (squeeze(CZ(1,1,:)) + squeeze(CZ(2,2,:)).*abs(squeeze(G(1,1,:).^2)) - ...
                      2*real(squeeze(CZ(1,2,:)).* conj(squeeze(G(1,1,:))))) .* (squeeze(SUInv(1,1,:)));
    
else                % MIMO case
    
    for kk = 1:F        
        dummy = [eye(ny), -squeeze(G(:,:,kk))] * squeeze(CZ(:,:,kk)) * [eye(ny), -squeeze(G(:,:,kk))]';
        CvecG(:, :, kk) = kron(squeeze(SUInv(:,:,kk)), dummy);        
    end % kk, frequencies
    
end % if one input
