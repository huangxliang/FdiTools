function SUinv = InvMeanInput(Zall)
%
%   function SUinv = InvMeanInput(Zall)
%
%
%   Output
%
%       SUinv   =   inverse of the mean (over the realisations) input matrix inv(conj(U*U')) 
%                   size nu x nu x F
%
%
%   Input
%
%       Zall    =   mean (over the periods) input DFT spectra  for all realisations
%                   size nu x nu x M x F
%
%   Rik Pintelon, October 29, 2009
%

[xx, nu, M, F] = size(Zall);

SUinv = zeros(nu, nu, F);

if nu == 1
    
    SUinv(1,1,:) = 1./squeeze(mean(abs(Zall.^2),3));
    
else % nu > 1
    
    % matrix product for all frequencies and realisations
    for kk = 1:F
        for mm = 1:M
            Zall(:,:,mm,kk) = conj(squeeze(Zall(:,:,mm,kk)) * squeeze(Zall(:,:,mm,kk))');
        end % mm, realisations
    end % kk, frequencies

    % mean value over the realisations + removal singular realisation dimension
    SUinv(:,:,:) = squeeze(mean(Zall, 3));

    % matrix inverse for all frequencies
    for kk = 1:F
        SUinv(:,:,kk) = inv(SUinv(:,:,kk));
    end % kk, frequencies

end % if one input
