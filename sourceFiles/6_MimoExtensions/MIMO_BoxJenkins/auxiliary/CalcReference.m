function data = CalcReference(data);
%
%   function data = CalcReference(data);
%
%
%	Output parameters
%		data				=	structure containing the non-parametric data required for the identification (see Input parameters)
%                               with additonal variable
%                                   data.R     =   DTF spectrum nu x 1 reference signal, dimensions: nu x number of frequencies
%
%
%	Input parameters
%		data				=	structure containing the non-parametric data required for the identification
%									data.Y		=	DFT spectrum ny x 1 output signal, dimensions ny x number of frequencies
%									data.U		=	DFT spectrum nu x 1 input signal, dimensions: nu x number of frequencies
%									data.freq	=	vector of frequency values (Hz), dimension: number of frequencies x 1
%									data.Ts		=	sampling time (s)
%									data.Gc		=	controller transfer function, dimension nu x ny x number of frequencies
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, April 2005 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%


%%%%%%%%%%%%%%%%%%%%%%%%%
% fast calculation Gc*Y %
%%%%%%%%%%%%%%%%%%%%%%%%%

[nu, ny, F] = size(data.Gc);
YT = zeros(1, ny, F);
YT(1,:,:) = data.Y;
GcY = sum(data.Gc.*repmat(YT, nu, 1), 2);
[rows, columns, depth] = size(GcY);
GcY = squeeze(GcY);
% squeeze function on 1 x 1 x F delivers F x 1 !!!
if rows*columns == 1
	GcY = GcY.'; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation R = U + Gc*Y %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.R = data.U + GcY;