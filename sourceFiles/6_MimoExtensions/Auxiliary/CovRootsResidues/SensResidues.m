function [S_Res, TheRoots] = SensResidues(B, A, ThePlane);
%
%   [S_Res, TheRoots] = SensResidues(B, A, ThePlane);
%
%
%               Calculates the sensitivity of the residue matrices of an ny x nu transfer function matrix G = B/A
%               w.r.t. to the poles of G, and the numerator coefficients of G;
%                   A is a polynomial of order na in increasing power of x
%                   B is an ny x nu matrix polynomial of order nb in increasing power of x
%                       s-, sqrt(s)-domain  :   x = s, sqrt(s)
%                       z-domain            :   x = z^-1
%
%
%   OUTPUT PARAMETER
%
%       S_Res           =   struct{'poles','B'} containing the sensitivity of the residues of the transfer function matrix G = B/A
%                           w.r.t. the poles of G, and the numerator coefficients of G contained in B
%                               S_Res.poles =   derivative residues w.r.t. ALL the poles of G; size ny x nu x na x na
%                               S_Res.B     =   derivative residues w.r.t. ALL the numerator coefficients of G; size ny x nu x ny x nu x (nb+1) x na
%                           Note: the residues are sorted in the same order as TheRoots
%
%       TheRoots        =   roots of the polynomial; size (number of roots) x 1
%                           Note: the roots are sorted in increasing order of magnitude
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
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny, nu, nb] = size(B);
nb = nb-1;
na = length(A)-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the residues and the poles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TheResidues, ThePoles] = CalcResidues(B, A, ThePlane);     % already sorted in increasing order of magnitude of ThePoles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impose a monic denominator polynomial %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ThePlane
    case {'s','w'}       
        for ii = 1:nb+1
            B(:,:,ii) = B(:,:,ii) / A(end);
        end % ii
        A = A / A(end);
    case 'z'    
        for ii = 1:nb+1
            B(:,:,ii) = B(:,:,ii) / A(1);
        end % ii
        A = A / A(1);
end % switch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity residues w.r.t. the poles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_Res.poles = zeros(ny, nu, na, na);

% derivative Res(ii) w.r.t. pole(jj) with ii <> jj
% S_Res.poles(:, :, jj, ii) = TheResidues(:, :, ii) / (ThePoles(ii) - ThePoles(jj));
if na > 1
	DiffPoles = zeros(1, 1, na);
	for jj = 1:na
        DiffPoles(1, 1, :) = ThePoles - ThePoles(jj);
        TheIndex = [1:jj-1, jj+1:na];
        S_Res.poles(:, :, jj, TheIndex) = TheResidues(:, :, TheIndex) ./ repmat(DiffPoles(TheIndex),[ny, nu, 1]);   
	end % jj
end % if

% derivative Res(ii) w.r.t. pole(ii)

	% calculation of the first term of the derivative of Res(ii) w.r.t. pole(ii)
	% new denominator coefficient: s-, sqrt(s)-domain: m*Bm; z-domain: (1-m)*Bm
	Bnew = zeros(ny, nu, nb+1);
 	TheIndex = zeros(1, 1, nb+1);
    switch ThePlane
        case {'s','w'}
            TheIndex(1, 1, :) = [0:nb];
        case 'z'
            TheIndex(1, 1, :) = 1 - [0:nb];
    end
    Bnew = B .* repmat(TheIndex, [ny, nu, 1]);
	FirstTerm = CalcResidues(Bnew, A, ThePlane);
	% required residue in the derivative is TheNewResidues/ThePoles
	dummy(1, 1, :) = ThePoles.';
	FirstTerm = FirstTerm ./ repmat(dummy, [ny, nu, 1]);

	% calculation of the second term of the derivative of Res(ii) w.r.t. pole(ii)
	switch ThePlane
        case {'s','w'}
            NewS_Res_poles = S_Res.poles;
        case 'z'
            TheMat = zeros(1, 1, na, na);
            TheMat(1, 1, :, :) = kron(ThePoles, (1./ThePoles).');
            NewS_Res_poles = S_Res.poles .* repmat(TheMat, [ny, nu, 1, 1]);
	end % switch
	SecondTerm = sum(permute(NewS_Res_poles, [1, 2, 4, 3]), 4);   % sum over the poles; permute is necessary otherwise a singleton dim remains

for ii = 1:na
    S_Res.poles(:, :, ii, ii) = FirstTerm(:, :, ii) - SecondTerm(:, :, ii);
end % ii


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity residues w.r.t. the numerator coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_Res.B = zeros(ny, nu, ny, nu, nb+1, na);

DenSens = zeros(1, na);
for ii = 1:na    
    switch ThePlane
        case {'s','w'}
            dummy = ThePoles(ii) - ThePoles;
        case 'z'
            dummy = 1 - ThePoles/ThePoles(ii);        
    end % switch
    DenSens(1, ii) = prod(dummy([1:ii-1,ii+1:na]));
end % ii

% powers of the poles
TheMat = zeros(nb+1, na);
switch ThePlane
    case {'s','w'}
        dummy = ThePoles.';
    case 'z'
        dummy = 1./ThePoles.';
end %switch
for ii = 1:nb+1
    TheMat(ii, :) = dummy.^(ii-1);
end % ii

% matrix defining the sensitivity
TheMat = TheMat ./ repmat(DenSens, [nb+1, 1]);
if ThePlane == 'z'
    TheMat = TheMat .* repmat(ThePoles.', [nb+1, 1]);
end % if

% sensitivity w.r.t. numerator coefficients
for ii = 1:ny
    for jj = 1:nu
        S_Res.B(ii, jj, ii, jj, :, :) = TheMat;
    end % jj
end % ii

