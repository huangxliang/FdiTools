function CS = EstimStochNLdist(CNL, Cn, M);
%
%       Calculates the covariance of the stochastic nonlinear distortions
%       as the difference between the total covariance and the noise
%       covariance. If the diaginal elements of CNL are smaller than that
%       of Cn then the corresponding (co-)variances of CS are set to zero.
%
%   function CS = EstimStochNLdist(CNL, Cn, M);
%
%
%   Output parameters
%
%       CS      =   covariance of the stochastic nonlinear distortions w.r.t. one MIMO experiment 
%
%
%   Input parameters
%
%       CNL     =   total covariance matrix
%
%       Cn      =   noise covariance matrix
%
%       M       =   number of random phase realisations used to calculate CNL and Cn
%
%
%   Rik Pintelon, November 4, 2009
%

TheSize = size(CNL);
Ndim = length(TheSize);                         % number of dimensions of the array 

%%%%%%%%%%%%%%%
% Estimate CS %
%%%%%%%%%%%%%%%

CS = M*(CNL - Cn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if the diagonal elements of CNL are smaller than those of Cn %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select the dimensions in excess of 2
if Ndim > 2
    ForLoopDim = prod(TheSize(3:end));
else
    ForLoopDim = 1;
end % if more than 2 dimensions

for kk = 1:ForLoopDim
    IndexNegative = find(abs(diag(squeeze(CNL(:,:,kk))))-abs(diag(squeeze(Cn(:,:,kk)))) < 0);
    NumberNegative = length(IndexNegative);
    for ii = 1:NumberNegative
        CS(IndexNegative(ii),:,kk) = 0;        % set the row no. IndexNegative(ii) to zero
        CS(:,IndexNegative(ii),kk) = 0;        % set the column no. IndexNegative(ii) to zero
    end % ii, negative diagonal elements CS
end % kk, linear index over all dimensions in excess of 2

