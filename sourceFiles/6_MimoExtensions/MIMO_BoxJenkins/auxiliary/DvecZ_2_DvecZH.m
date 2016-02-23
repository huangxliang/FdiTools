function DvecZH = DvecZ_2_DvecZH(DvecZ, rows, columns);
%
% function DvecZH = DvecZ_2_DvecZH(DvecZ, rows, columns);
%
%       Given the derivative of vec(Z) w.r.t. to a vector param, calculate the derivative of vec(Z^H) w.r.t. param
%       where Z has dimensions rows x columns x number of frequencies.
%       Note: size(DvecZ, 1) = rows*columns
%
%   Output parameter
%       DvecZH  =   derivative of vec(Z') w.r.t. parameter vector, dimension rows*columns x dim(param) x number of frequencies
%                   where Z has dimensions rows x columns x number of frequencies
%
%   Input parameter
%       DvecZ   =   derivative of vec(Z) w.r.t. parameter vector, dimension rows*columns x dim(param) x number of frequencies
%                   where Z has dimensions rows x columns x number of frequencies
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, April 2005 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%

F = size(DvecZ, 3);
np = size(DvecZ, 2);
DvecZH = zeros(size(DvecZ));

% intermediate matrix with singleton dimension to calculate DvecZH for each column of DvecZ seperately
dummy = zeros(rows*columns, 1, np, F); 
dummy(:,1,:,:) = DvecZ;                          % derivative vec(Z) w.r.t. param
dummy = reshape(dummy, rows, columns, np, F);    % derivative Z w.r.t. param
dummy = conj(permute(dummy,[2, 1, 3, 4]));       % derivative Z' = Z^H w.r.t. param
dummy = reshape(dummy, rows*columns, 1, np, F);  % derivative vec(Z^H) w.r.t. param
dummy = permute(dummy, [1, 3, 4, 2]);

DvecZH = dummy;     
