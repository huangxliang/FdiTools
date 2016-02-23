function CRbound = CRtransferfunctions(CRbound, SqrtCRtheta, PolyTrans, Deriv, Sel, ModelVar);
%
% function CRbound = CRtransferfunctions(CRbound, SqrtCRtheta, PolyTrans, Deriv, Sel, ModelVar);
%
%
%   Output parameter
%
%		CRbound		=	see input parameter; the following fields are added 'vecG' and 'vecH' which are the
%                           Cramer-Rao lower bounds of the plant (vec(G)) and noise (vec(H)) transfer functions
%							CRbound.vecG = ny*nu x ny*nu x number of freq.
%								CRbound.vecG(i,j) = covariance between plant transfer functions vecG(i) and vecG(j) 
%							CRbound.G = ny x nu x F
%								CRbound.G(i,j,r) = variance G(i,j,r)
%							CRbound.vecH = (ny^2) x (ny^2) x F
%								CRbound.vecH(i,j,r) = covariance between vecH(i,r) and vecH(j,r)
%							CRbound.H = ny x ny x F
%								CRbound.H(i,j,r) = variance H(i,j,r)
%
%   Input parameters
%
%		CRbound		=	Cramer-Rao bound of the estimated model parameters, the estimated plant model, and the estimated noise model
%						structure with fields 'A', 'vecB', 'vecC', 'D', 'Theta'
%							CRbound = struct('A', [], 'vecB',[], 'vecC', [], 'D', [], 'Theta', [])
%							CRbound.A = FreeParam.A x FreeParam.A
%								CRbound.A(i,j) = covariance between free coefficients a(i-1) and a(j-1) 
%							CRbound.vecB = FreeParam.B x FreeParam.B
%								CRbound.vecB(i,j) = covariance between free parameters vecB(i) and vecB(j)
%							CRbound.vecC = FreeParam.C x FreeParam.C
%								CRbound.vecC(i,j) = covariance between free parameters vecC(i) and vecC(j) 
%							CRbound.D = FreeParam.D x FreeParam.D
%								CRbound.D(i,j) = covariance between free coefficients d(i-1) and d(j-1) 
%							CRbound.all = dim(Theta) x dim(Theta), where dim(Theta) = FreeParam.A + FreeParam.B + FreeParam.C + FreeParam.D
%								CRbound.all(i,j) = covariance between free parameters Theta(i) and Theta(j)
%							Notes:
%                               - in s-, sqrt(s) domains the CR-bound of the normalised parameters is calculated
%                               - the (normalised) model parameters satisfy the following constraints
%                                   z-domain:               a(0) = d(0) = 1
%                                   s-, sqrt(s)-domains:    a(OrderA) = d(OrderD) = 1
%                               - the C-coefficients are normalised such that CovE = eye(ny)
%       SqrtCRtheta =   square root of CRbound.all to guarantee a numerical stable calculation of CRbound.vecG, CRbound.vecH, ...
%                       CRbound.all = SqrtCRtheta * SqrtCRtheta.'
%		PolyTrans	=	structure containing the polynomials and transfer functions evaluated in x
%							PolyTrans.A		=	denominator polynomial plant transfer function evaluated in x.Plant, dimensions 1 x number of freq.
%							PolyTrans.D		=	D polynomial evaluated in x.Noise, dimensions 1 x number of freq.
%							PolyTrans.G		=	plant transfer matrix evaluated in x.Plant, dimensions ny x nu x number of freq.
%							PolyTrans.H		=	noise transfer matrix evaluated in x.Noise, dimensions ny x ny x number of freq.
%							PolyTrans.Hinv	=	inverse of the noise transfer matrix evaluated in x.Noise, dimensions ny x ny x number of freq.
%							PolyTrans.Tg	=	plant transient term evaluated in x.Plant, dimension ny x number of freq.
%							PolyTrans.Th	=	noise transient term evaluated in x.Plant, dimension ny x number of freq.
%		Deriv		=	structure containing the derivative of vec(G), vec(H), vec(G'), vec(H') w.r.t.
%						all the plant model parameters a, b, c, d
%						    Deriv.vecGa	    =	derivative vec(G) w.r.t. a;	size ny*nu x (na+1) x number of freq.
%						    Deriv.vecGb	    =	derivative vec(G) w.r.t. b;	size ny*nu x ny*nu*(nb+1) x number of freq.
%						    Deriv.vecHc	    =	derivative vec(H) w.r.t. c;	size ny^2 x ny^2*(nc+1) x number of freq.
%						    Deriv.vecHd	    =	derivative vec(H) w.r.t. d;	size ny^2 x (nd+1) x number of freq.
%						    Deriv.vecGHa	=	derivative vec(G') w.r.t. a; size ny*nu x (na+1) x number of freq.
%						    Deriv.vecGHb	=	derivative vec(G') w.r.t. b; size ny*nu x ny*nu*(nb+1) x number of freq.
%						    Deriv.vecHHc	=	derivative vec(H') w.r.t. c; size ny^2 x ny^2*(nc+1) x number of freq.
%						    Deriv.vecHHd	=	derivative vec(H') w.r.t. d; size ny^2 x (nd+1) x number of freq.
%		Sel			=	structure with fields 'A', 'B', 'Ig', 'C', 'D', 'Ih'
%							Sel = struct('A',[],'B',[], 'Ig', [], 'C',[],'D',[], 'Ih', [])
%							Sel.A = 1 x (OrderA+1)
%								Sel.A(r) = 1 if coeff. a(r-1) is unknown
%								Sel.A(r) = 0 if coeff. a(r-1) = 0
%							Sel.B = ny x nu x (OrderB+1)
%								Sel.B(i,j,r) = 1 if coeff. b(i,j,r-1) is unknown
%								Sel.B(i,j,r) = 0 if coeff. b(i,j,r-1) = 0
%							Sel.Ig = ny x (OrderIg+1)
%								Sel.Ig(i,r) = 1 if coeff. ig(i,r-1) is unknown
%								Sel.Ig(i,r) = 0 if coeff. ig(i,r-1) = 0
%							Sel.C = ny x ny x (OrderC+1)
%								Sel.C(i,j,r) = 1 if coeff. c(i,j,r-1) is unknown
%								Sel.C(i,j,r) = 0 if coeff. c(i,j,r-1) = 0
%							Sel.D = 1 x (OrderD+1)
%								Sel.D(r) = 1 if coeff. d(i,j,r-1) is unknown
%								Sel.D(r) = 0 if coeff. d(i,j,r-1) = 0
%							Sel.Ih = ny x (OrderIh+1)
%								Sel.Ih(i,r) = 1 if coeff. ih(i,r-1) is unknown
%								Sel.Ih(i,r) = 0 if coeff. ih(i,r-1) = 0
%		ModelVar	=	contains the information about the model to be identified structure with the following fields
%							ModelVar.Transient		=	1 then the initial conditions of the plant and/or noise are estimated
%							ModelVar.PlantPlane		=	plane of the plant model
%															's':	    continuous-time
%															'w':	    sqrt(s)-domain
%															'z':	    discrete-time
%															'':		    plane not defined
%							ModelVar.NoisePlane		=	plane of the plant model
%															's':	    continuous-time
%															'w':	    sqrt(s)-domain
%															'z':	    discrete-time
%															'':		    plane not defined
%							ModelVar.Struct			=	model structure
%															'BJ':		Box-Jenkins
%															'OE':		output error (plant model only)
%															'ARMA':		autoregressive moving average (noise model only)
%															'ARMAX':    autoregressive moving average with exogenous input
%							ModelVar.RecipPlant		=	1 if plant model is reciprocal: G(i,j) = G(j,i)
%							ModelVar.nu				=	number of inputs
%							ModelVar.ny				= 	number of outputs
%							ModelVar.na				=	order polynomial A
%							ModelVar.nb				= 	order matrix polynomial B
%							ModelVar.nc				=	order matrix polynomial C
%							ModelVar.nd				=	order polynomial D
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, April 2005 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version June, 2005
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix model parameters:                                                %
%   1. collect the derivatives of vec(G) and vec(H) w.r.t. theta into one matrix     %
%   2. put the derivative matrices in the following form:                            %
%       nu*ny x ntheta x F => nu*ny*F x ntheta for the derivatives of vec(G)         %
%       ny^2 x ntheta x F => ny^2*F x ntheta for the derivatives of vec(H)           %
%   3. impose the common parameter structure and eliminate excess parameters         %
%   4. put the derivative matrices back in their original form:                      %
%       nu*ny*F x nfreetheta => nu*ny x nfreetheta x F for the derivatives of vec(G) %
%       ny^2*F x nfreetheta => ny^2 x nfreetheta x F for the derivatives of vec(H)   %
%   5. calculate the covariance matrices for all frequencies                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu = ModelVar.nu;
ny = ModelVar.ny;
na = ModelVar.na;
nb = ModelVar.nb;
nig = ModelVar.nig;
nc = ModelVar.nc;
nd = ModelVar.nd;
nih = ModelVar.nih;
ntheta = (na+1) + (nb+1)*nu*ny + (nig+1)*ny + (nc+1)*ny^2 + (nd+1) + (nih+1)*ny; % total number of model parameters
F = size(PolyTrans.G, 3);

% 1. collect the derivatives of vec(G) and vec(H) w.r.t. ALL model parameters (without constraints) into one matrix
DerivVecG = [Deriv.vecGa, Deriv.vecGb, zeros(nu*ny, (nig+1)*ny + (nc+1)*ny^2 + (nd+1) + (nih+1)*ny, F)];
DerivVecH = [zeros(ny^2, (na+1) + (nb+1)*nu*ny + (nig+1)*ny, F), Deriv.vecHc, Deriv.vecHd, zeros(ny^2, (nih+1)*ny, F)];

% 2. reshape the derivative matrices
%    the number of columns of the derivative matrices equals the number
%    of ALL model parameters (without constraints)
DerivVecG = reshape(permute(DerivVecG, [1, 3, 2]), [nu*ny*F, ntheta]);
DerivVecH = reshape(permute(DerivVecH, [1, 3, 2]), [ny^2*F, ntheta]);

% 3. impose the common parameter structure and eliminate the excess parameters
%    the number of columns of the derivative matrices equals then the number
%    of FREE model parameters
DerivVecG = Add_SelectColumns(DerivVecG, Sel, ModelVar);
DerivVecH = Add_SelectColumns(DerivVecH, Sel, ModelVar);

% 4. reshape the derivative matrices to their original form
FreeParam.Theta = size(DerivVecG, 2);
DerivVecG = permute(reshape(DerivVecG, [nu*ny, F, FreeParam.Theta]), [1, 3, 2]);
DerivVecH = permute(reshape(DerivVecH, [ny^2, F, FreeParam.Theta]), [1, 3, 2]);

% 5. calculation of the covariance matrices of the plant and noise transfer functions
CRbound.vecG = zeros(nu*ny, nu*ny, F);
CRbound.G = zeros(ny, nu, F);
CRbound.vecH = zeros(ny^2, ny^2, F);
CRbound.H = zeros(ny, ny, F);
for kk = 1:F
    XX = DerivVecG(:, :, kk) * SqrtCRtheta;
    CRbound.vecG(:, :, kk) = XX * XX';
    CRbound.G(:, :, kk) = reshape(diag(CRbound.vecG(:, :, kk)), [ny, nu]);
    XX = DerivVecH(:, :, kk) * SqrtCRtheta;
    CRbound.vecH(:, :, kk) = XX * XX';
    CRbound.H(:, :, kk) = reshape(diag(CRbound.vecH(:, :, kk)), [ny, ny]);
%     CRbound.vecG(:, :, kk) = DerivVecG(:, :, kk) * CRbound.all * DerivVecG(:, :, kk)';
%     CRbound.vecG(:, :, kk) = CRbound.vecG(:, :, kk) - sqrt(-1) * diag(imag(diag(CRbound.vecG(:, :, kk))));
%     CRbound.G(:, :, kk) = reshape(diag(CRbound.vecG(:, :, kk)), [ny, nu]);
%     CRbound.vecH(:, :, kk) = DerivVecH(:, :, kk) * CRbound.all * DerivVecH(:, :, kk)';
%     CRbound.vecH(:, :, kk) = CRbound.vecH(:, :, kk) - sqrt(-1) * diag(imag(diag(CRbound.vecH(:, :, kk))));
%     CRbound.H(:, :, kk) = real(reshape(diag(CRbound.vecH(:, :, kk)), [ny, ny]));
end % kk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix noise power spectra and their square root %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% derivative vec(H*H') w.r.t. the free model parameters
DerivVecHH = DvecZ_2_DvecZH(DerivVecH, ny, ny);
DerivNoisePower = zeros(ny^2, FreeParam.Theta, F);
for kk = 1:F
    DerivNoisePower(:, :, kk) = kron(conj(PolyTrans.H(:, :, kk)), eye(ny)) * DerivVecH(:, :, kk) + ...
                                kron(eye(ny), PolyTrans.H(:, :, kk)) * DerivVecHH(:, :, kk);
end % kk

% Cramer-Rao bound entries noise power spectrum
% in an intermediate step var(vec(H*H')), size ny^2 x F
CRbound.NoisePower = zeros(ny^2, F);
NoisePower = zeros(ny, ny, F);
for kk = 1:F
    CRbound.NoisePower(:, kk) = diag(DerivNoisePower(:, :, kk) * CRbound.all * DerivNoisePower(:, :, kk)');
    NoisePower(:, :, kk) = PolyTrans.H(:, :, kk) * PolyTrans.H(:, :, kk)';
end %kk
% reshape to var(H*H'), size ny x ny x F
CRbound.NoisePower = reshape(CRbound.NoisePower, [ny, ny, F]);
if ~strcmp(ModelVar.Struct,'OE')
    CRbound.SqrtNoisePower = CRbound.NoisePower./(4*abs(NoisePower));
else
    CRbound.SqrtNoisePower = [];
end
