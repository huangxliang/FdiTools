function [y, z] = WH_LPM(u, Gfirst, fNL, scale, Gsecond)
%
%       Calculates the steady state response of a continuous-time
%       Wiener-Hammerstein system to a periodic input 
%
% 		[y, z] = WH_LPM(u, Gfirst, fNL, scale, Gsecond, u);
%
%
%	output parameters
%
%		y		: the time signal at the output of the Wiener-Hammerstein block
%		z		: the time signal at the output of the static nonlinearity (to verify if the
%				  the oversampling is large enough
%
%	input parameters
%
%		u		: one period of the input time signal
%		Gfirst	: transfer function of the first linear dynamic block from DC to Nyquist.
%		fNL		: string containing the (m-file) function describing the static nonlinear block
%		scale	: input fNL is scaled by scale, and output fNL is multipied by scale 
%		Gsecond	: transfer function of the second linear dynamic block from DC to Nyquist.
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 3 December 2009 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 3 October 2011
%

% number of samples
N = length(u);

% signal after first linear dynamic system: X = Gfirst*U
U = fft(u);
X = zeros(1,N);
X(1:N/2+1) = Gfirst.*U(1:N/2+1);
% take 2 times real part since complex conjugate was not added in X
% Note: DC and Nyquist should not be scaled by 2
X(2:N/2) = 2*X(2:N/2);
x = real(ifft(X));

% signal after static nonlinearity: z = scale*fNL(x/scale)
z = scale*feval(fNL,x/scale);

% signal after second linear dynamic system: Y = Gsecond*Z
Z = fft(z);
Y = zeros(1,N);
Y(1:N/2+1) = Gsecond.*Z(1:N/2+1);
% take 2 times real part since complex conjugate was not added in Y
% Note: DC and Nyquist should not be scaled by 2
Y(2:N/2) = 2*Y(2:N/2);
y = real(ifft(Y));