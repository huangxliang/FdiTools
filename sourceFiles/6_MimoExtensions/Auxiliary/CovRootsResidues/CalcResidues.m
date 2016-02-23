function [TheResidues, ThePoles] = CalcResidues(B, A, ThePlane);
%
%        [TheResidues, ThePoles] = CalcResidues(B, A, ThePlane);
%
%               Calculates the residue matrices and the poles of a ny x nu transfer function matrix G = B/A
%               with A a polynomial of order na, and B an ny x nu matrix polynomial of order nb
%
%   OUTPUT PARAMETERS
%
%       TheResidues     =   residue matrices of the ny x nu transfer function matrix G = B/A; size ny x nu x na
%                           Note: the residues are in the same order of ThePoles
%
%       ThePoles        =   poles of the transfer function matrix; size na x 1
%                           Note: the poles are sorted in increasing order of magnitude
%
%
%   INPUT PARAMETERS
%
%       B               =   matrix polynomial of order nb in increasing powers of s, sqrt(s), or z^-1; size ny x nu x (nb+1)
%
%       A               =   denominator polynomial of order na in increasing powers of s, sqrt(s), or z^-1; size 1 x (na+1)
%
%       ThePlane        =   domain of the transfer function matrix
%                               z-domain:       'z'
%                               s-domain:       's'
%                               sqrt(s)-domain: 'w'
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, February 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 6 October 2011
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny, nu, nb] = size(B);
nb = nb-1;
na = length(A)-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the denominator polynomial monic and calculate the poles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ThePlane
    case {'s','w'}       
        for ii = 1:nb+1
            B(:,:,ii) = B(:,:,ii)/A(end);
        end % ii
        A = A/A(end);
        ThePoles = sort(roots(fliplr(A)));
    case 'z'    
        for ii = 1:nb+1
            B(:,:,ii) = B(:,:,ii)/A(1);
        end % ii
        A = A/A(1);
        ThePoles = sort(roots(A));
end % switch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the residue matrices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% numerator residue matrices
NumRes = zeros(ny, nu, na);
for ii = 1:ny 
    for jj = 1:nu
        switch ThePlane
            case {'s','w'}
                NumRes(ii, jj, :) = polyval(fliplr(squeeze(B(ii, jj, :)).'), ThePoles.');
            case 'z'
                NumRes(ii, jj, :) = polyval(fliplr(squeeze(B(ii, jj, :)).'), (1./ThePoles).');
        end % switch
    end % jj
end % ii


% denominator residue matrices
DenRes = zeros(1, na);
for ii = 1:na    
    switch ThePlane
        case {'s','w'}
            dummy = ThePoles(ii) - ThePoles;
        case 'z'
            dummy = 1 - ThePoles/ThePoles(ii);        
    end % switch
    DenRes(1, ii) = prod(dummy([1:ii-1,ii+1:na]));
end % ii
if ThePlane == 'z'
    DenRes = DenRes ./ (ThePoles.');
end


% residue matrices
dummy = zeros(1, 1, na);
dummy(1, 1, :) = DenRes;
TheResidues = NumRes ./ repmat(dummy, [ny, nu, 1]);

