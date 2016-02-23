function [Gc, Grat] = CorrectTimeVariationFRF(G);
%
%       function Gc = CorrectTimeVariationFRF(G);
%
%       First order correction of the time variation of the FRF
%
%
%   OUTPUT
%
%       Gc      =   structure {'all', 'mean', 'stdn', 'stdNL'} containing mean value and standard deviations of the FRF at excited odd harmonics,
%                   with a first order correction for the time variation
%                   Gc.all      =   Frequency Response Function (FRF) for all realisations and periods; size M x P x F
%                   Gc.mean     =   mean value FRF over the P consecutive periods
%                   Gc.stdn     =   struct{'E', 'NE'} containing the noise standard deviation of the mean FRF value over the P consecutive periods
%                                   Gc.stdn.E    =   noise std calculated from the excited frequencies
%                                   Gc.stdn.NE   =   noise std calculated from the non-excited odd frequencies (= extrapolation)
%                                   Note: a difference between G.stdn.E and G.stdn.NE indicates a non-stationary behaviour
%                   Gc.stdNL    =   total standard deviation mean value FRF over P consecutive periods:
%                                   Gc.stdNL.^2  =   noise variance + stochatic NL distortions; size M x F
%
%       Grat    =   struct{'all', 'mean'} containg the ratios of the FRF's over consecutive periods
%                   Grat.all    =   ratio Gall(:, ii+1, :) ./ Gall(:, ii, :)
%                   Grat.mean   =   mean value Grat over consecutive periods
%
%   INPUT
%
%       G       =   structure {'all', 'mean', 'stdn', 'stdNL'} containing mean value and standard deviations FRF at excited odd harmonics;
%                   G.all       =   Frequency Response Function (FRF) for all realisations and periods; size M x P x F or P x F
%                   G.mean      =   mean value FRF over the P consecutive periods; size M x F
%                   G.stdn      =   struct{'E', 'NE'} containing the noise standard deviation of the mean FRF value over the P consecutive periods
%                                   G.stdn.E    =   noise std calculated from the excited frequencies; size M x F
%                                   G.stdn.NE   =   noise std calculated from the non-excited odd frequencies (= extrapolation); size M x F
%                                   Note: a difference between G.stdn.E and G.stdn.NE indicates a non-stationary behaviour
%                   G.stdNL     =   total standard deviation mean value FRF over P consecutive periods:
%                                   G.stdNL.^2  =   noise variance + stochatic NL distortions; size M x F
%
%
%   Rik Pintelon, April 2006
%

% initialisation variables
TheSize = size(G.all);
if length(TheSize) == 3
    M = TheSize(1);
    P = TheSize(2);
    F = TheSize(3);
else
    M = 1;
    P = TheSize(1);
    F = TheSize(2);
    dummy = zeros(1, P, F);
    dummy(1, :, :) = G.all;
    G.all = dummy;
end % if
Gc = G;

% ratio FRF consecutive periods
Grat.all = G.all(:, 2:P, :) ./ G.all(:, 1:P-1, :);
Grat.mean = mean(Grat.all, 2);

% first order correction time variation FRF
TheScale = zeros(M, P, F);
TheScale(:, 1, :) = ones(M, F);
for ii = 1:P-1
    TheScale(:, ii+1, :) = Grat.mean.^ii;
end % ii
Gc.all = G.all ./ TheScale;
Gc.stdn.E = std(Gc.all, 0, 2) / sqrt(P);
if M == 1
    Gc.stdn.E = squeeze(Gc.stdn.E).';
end % if

