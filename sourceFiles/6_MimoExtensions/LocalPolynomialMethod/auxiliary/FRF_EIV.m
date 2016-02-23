function [Geiv, CvecGeiv] = FRF_EIV(G, CvecG);
%
%		Calculates the input/output FRF and its covariance matrix from the reference signal 
%       to the input/output FRF and its covariance matrix.
%       (Note: Reference signal: nu x 1; input signal: nu x 1; output signal: ny x 1) 
%       
%
%	function [Geiv, CvecGeiv] = FRF_EIV(G, CvecG);
%
%	Output parameters
%
%
%       Geiv        =   estimated input/output frequency response matrix, size ny x nu x F
%
%       CvecGeiv	=   covariance matrix vec(G), size (ny*nu) x (ny*nu) x F
%                       or struct('NL', [], 'n', []) containing the total and noise covariance of vec(Geiv) 
%                           CvecGeiv.NL     =   total covariance of vec(Geiv), size ((ny+nu)*nu) x ((ny+nu)*nu) x F 
%                           CvecGeiv.n      =   noise covariance of vec(Geiv), size ((ny+nu)*nu) x ((ny+nu)*nu) x F 
%
%	Input parameters
%
%       G           =   estimated frequency response matrix, size (ny + nu) x nu x F
%
%       CvecG       =   covariance matrix vec(G), size ((ny+nu)*nu) x ((ny+nu)*nu) x F 
%                       or struct('NL', [], 'n', []) containing the total and noise covariance of vec(G) 
%                           CvecG.NL        =   total covariance of vec(G), size ((ny+nu)*nu) x ((ny+nu)*nu) x F 
%                           CvecG.n         =   noise covariance of vec(G), size ((ny+nu)*nu) x ((ny+nu)*nu) x F 
%
%
% Rik Pintelon, November 21, 2008
% version January 4, 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    if isstruct(CvecG)
        StructCov = 1;          % CvecG is a structure of matrices
    else
        StructCov = 0;          % CvecG is a matrix
    end
catch
    StructCov = 0;              % CvecG is a matrix
end % try

[nz, nu, F] = size(G);
ny = nz-nu;                     % number of outputs

Geiv = zeros(ny, nu, F);
if StructCov == 0               % CvecG is a matrix
    CvecGeiv = zeros(ny*nu, ny*nu, F);
else                            % CvecG is a structure of matrices
    CvecGeiv = struct('NL', [], 'n', []);
    CvecGeiv.NL = zeros(ny*nu, ny*nu, F);
    CvecGeiv.n = zeros(ny*nu, ny*nu, F);    
end % if CvecG is a matrix


%%%%%%%%%%%%%%%%%%%%
% Input-output FRF %
%%%%%%%%%%%%%%%%%%%%

for kk = 1:F
    
    % inverse FRF from reference to input
    Gru_inv = inv(squeeze(G(ny+1:end,:,kk)));
 
    % input-output FRF
    Geiv(:,:,kk) = G(1:ny,:,kk) * Gru_inv;
    
    % covariance matrix input-output FRF
    dummy = kron(Gru_inv.', [eye(ny), -squeeze(Geiv(:,:,kk))]); 
    if StructCov == 0               % CvecG is a matrix
        CvecGeiv(:,:,kk) = dummy * CvecG(:,:,kk) * dummy';
    else                            % CvecG is a structure of matrices
        CvecGeiv.NL(:,:,kk) = dummy * CvecG.NL(:,:,kk) * dummy';
        CvecGeiv.n(:,:,kk) = dummy * CvecG.n(:,:,kk) * dummy';
    end % if CvecG is a matrix
    
end % for kk
