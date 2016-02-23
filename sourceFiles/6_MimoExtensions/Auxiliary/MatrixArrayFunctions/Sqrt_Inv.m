function y = Sqrt_Inv(x, delta);
%
%   function y = Sqrt_Inv(x, delta);
%
%
%   Output
%
%       y   =   symmetric square root of the inverse of the first two matrix dimensions of x 
%
%
%   Input
%
%       x   =   nx x nx x F matrix array of hermitian symmetric positive definite matrices 
%
%       delta   =   relative lower limit under which the singular values are considered to be zero 
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 1 December 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%

[nx, nx, F] = size(x);
y = zeros(nx,nx,F);

for kk = 1:F
    
    [ux, sx, vx] = svd(squeeze(x(:,:,kk)),0);
    sx = diag(sx).^0.5;
    Index = find(sx/sx(1) <= delta);
    sx(Index) = inf;
    y(:,:,kk) = ux * diag(1./sx) * ux';
    
end % kk, frequencies