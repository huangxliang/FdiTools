function [S_RootsCoeff, TheRoots] = SensRootsCoeff(Coeff, ThePlane);
%
%   [S_RootsCoeff, TheRoots] = SensRootsCoeff(Coeff, ThePlane);
%
%
%       Calculates the sensitivity matrix of the roots of a polynomial w.r.t. its coefficients
%
%               Coeff(1) + Coeff(2)*x + Coeff(3) * x^2 + ... + Coeff(nx+1) x^nx
%
%       s-, sqrt(s)-domain  :   x = s, sqrt(s)
%       z-domain            :   x = z^-1
%
%
%   OUTPUT PARAMETER
%
%       S_RootsCoeff    =   sensitivity of the roots w.r.t. the polynomial coefficients; size (number of roots) x number of coefficients
%
%       TheRoots        =   roots of the polynomial; size (number of roots) x 1
%                           Note: the roots are sorted in increasing order of magnitude
%
%
%   INPUT PARAMETERS
%
%       Coeff           =   coefficients of a monic polynomial in increasing powers of s, sqrt(s), or z^-1; size 1 x (1 + order polynomial)
%
%       ThePlane        =   domain of the polynomial
%						        's':	continuous-time;
%						        'w':	sqrt(s)-domain
%						        'z':	discrete-time;
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, January 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%


%%%%%%%%%%%%%%%%%%
% initialisation %
%%%%%%%%%%%%%%%%%%

% order polynomial
nn = length(Coeff) - 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity of the roots w.r.t. the estimated polynomial coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the roots
% in the z-domain 1/TheRoots is calculated
TheRoots = roots(fliplr(Coeff));

% coefficient of the differentiated polynomial
dCoeff = [1:nn];
dCoeff = Coeff(2:end) .* dCoeff;

% differentiated polynomial evaluated in the roots of the original polynomial
dPoly = polyval(fliplr(dCoeff), TheRoots);

% table containing the powers of the roots
PowerRoots = zeros(nn, nn+1);
for ii = 1:nn+1
    PowerRoots(:, ii) = TheRoots.^(ii-1);    
end

% sensitivity matrix of the roots w.r.t. the estimated polynomial coefficients
switch ThePlane
    case {'s', 'w'}
        S_RootsCoeff = - PowerRoots ./ repmat(dPoly, 1, nn+1);
    case 'z'
        S_RootsCoeff = PowerRoots ./ repmat(dPoly .* TheRoots.^2, 1, nn+1);
end
 
% roots in z-domain
if ThePlane == 'z'
    TheRoots = 1./TheRoots;    
end

% sort the roots in increasing order of magnitude
[TheRoots, Index] = sort(TheRoots);
S_RootsCoeff = S_RootsCoeff(Index, :);