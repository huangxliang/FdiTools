function CRbound = MIMO_ML_CRtf(CRbound, SqrtCRtheta, PolyTrans, Deriv, Sel, ModelVar);
%
% function CRbound = MIMO_ML_CRtf(CRbound, SqrtCRtheta, PolyTrans, Deriv, Sel, ModelVar);
%
%
%   Output parameter
%
%		CRbound		=	see input parameter; the following fields are added 'vecG' and 'vecH' which are the
%                           Cramer-Rao lower bounds of the plant (vec(G)) transfer function matrix 
%							CRbound.vecG = ny*nu x ny*nu x F
%								CRbound.vecG(i,j) = covariance between plant transfer functions vecG(i) and vecG(j) 
%							CRbound.G = ny x nu x F
%								CRbound.G(i,j,r) = variance G(i,j,r)
%
%
%   Input parameters
%
%		CRbound		=	Cramer-Rao bound of the estimated model parameters, the estimated plant model, and the estimated noise model
%						structure with fields 'A', 'vecB', 'Theta'
%							CRbound = struct('A', [], 'vecB', [], 'Theta', [])
%							CRbound.A = FreeParam.A x FreeParam.A
%								CRbound.A(i,j) = covariance between free coefficients a(i-1) and a(j-1) 
%							CRbound.AvecB = FreeParam.A x FreeParam.B
%								CRbound.AvecB(i,j) = covariance between free coefficients a(i-1) and vecB(j) 
%							CRbound.vecB = FreeParam.B x FreeParam.B
%								CRbound.vecB(i,j) = covariance between free parameters vecB(i) and vecB(j)
%							CRbound.all = dim(Theta) x dim(Theta), where dim(Theta) = FreeParam.A + FreeParam.B + FreeParam.C + FreeParam.D
%								CRbound.all(i,j) = covariance between free parameters Theta(i) and Theta(j)
%							Notes:
%                               - in s-, sqrt(s) domains the CR-bound of the normalised parameters is calculated
%                               - the (normalised) model parameters satisfy the following constraints
%                                   z-domain:               a(0)      = 1
%                                   s-, sqrt(s)-domains:    a(OrderA)  = 1
%
%       SqrtCRtheta =   square root of CRbound.all to guarantee a numerical stable calculation of CRbound.vecG, ...
%                       CRbound.all = SqrtCRtheta * SqrtCRtheta.' 
%
%		PolyTrans	=	structure containing the polynomials and transfer functions evaluated in x
%							PolyTrans.A		=	denominator polynomial plant transfer function evaluated in x.Plant, size 1 x F 
%							PolyTrans.G		=	plant transfer matrix evaluated in x.Plant, size ny x nu x F 
%							PolyTrans.Tg	=	plant transient term evaluated in x.Plant, size ny x F 
%
%		Deriv		=	structure containing the derivative of vec(G), vec(G') w.r.t. all the plant model parameters a, b
%						    Deriv.vecGa	    =	derivative vec(G) w.r.t. a; size ny*nu x (na+1) x F 
%						    Deriv.vecGb	    =	derivative vec(G) w.r.t. b;	size ny*nu x ny*nu*(nb+1) x F
%						    Deriv.vecGHa	=	derivative vec(G') w.r.t. a; size ny*nu x (na+1) x F
%						    Deriv.vecGHb	=	derivative vec(G') w.r.t. b; size ny*nu x ny*nu*(nb+1) x F
%
%		Sel			=	structure with fields 'A', 'B', 'Ig'
%							Sel = struct('A',[],'B',[], 'Ig', [])
%							Sel.A = 1 x (OrderA+1)
%								Sel.A(r) = 1 if coeff. a(r-1) is unknown
%								Sel.A(r) = 0 if coeff. a(r-1) = 0
%							Sel.B = ny x nu x (OrderB+1)
%								Sel.B(i,j,r) = 1 if coeff. b(i,j,r-1) is unknown
%								Sel.B(i,j,r) = 0 if coeff. b(i,j,r-1) = 0
%							Sel.Ig = ny x (OrderIg+1)
%								Sel.Ig(i,r) = 1 if coeff. ig(i,r-1) is unknown
%								Sel.Ig(i,r) = 0 if coeff. ig(i,r-1) = 0
%
%		ModelVar	=	contains the information about the model to be identified structure with the following fields
%							ModelVar.Transient		=	1 then the initial conditions of the plant and/or noise are estimated
%							ModelVar.PlantPlane		=	plane of the plant model
%															's':	continuous-time
%															'w':	sqrt(s)-domain
%															'z':	discrete-time
%															'':		plane not defined
%							ModelVar.Struct			=	model structure
%                                                           'EIV':	errors-in-variables (noisy input-output data)
%                                                           'OE':	generalised output error (known input, noisy output)
%							ModelVar.RecipPlant		=	1 if plant model is reciprocal: G(i,j) = G(j,i)
%							ModelVar.nu				=	number of inputs
%							ModelVar.ny				= 	number of outputs
%							ModelVar.na				=	order polynomial A
%							ModelVar.nb				= 	order matrix polynomial B
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix model parameters:                                                %
%   1. collect the derivatives of vec(G) w.r.t. theta into one matrix                %
%   2. put the derivative matrices in the following form:                            %
%       nu*ny x ntheta x F => nu*ny*F x ntheta for the derivatives of vec(G)         %
%   3. impose the common parameter structure and eliminate excess parameters         %
%   4. put the derivative matrices back in their original form:                      %
%       nu*ny*F x nfreetheta => nu*ny x nfreetheta x F for the derivatives of vec(G) %
%   5. calculate the covariance matrices for all frequencies                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu = ModelVar.nu;
ny = ModelVar.ny;
na = ModelVar.na;
nb = ModelVar.nb;
nig = ModelVar.nig;
ntheta = (na+1) + (nb+1)*nu*ny + (nig+1)*ny; % total number of model parameters
F = size(PolyTrans.G, 3);

% 1. collect the derivatives of vec(G) w.r.t. ALL model parameters (without constraints) into one matrix
DerivVecG = [Deriv.vecGa, Deriv.vecGb, zeros(nu*ny, (nig+1)*ny, F)];

% 2. reshape the derivative matrices
%    the number of columns of the derivative matrices equals the number
%    of ALL model parameters (without constraints)
DerivVecG = reshape(permute(DerivVecG, [1, 3, 2]), [nu*ny*F, ntheta]);

% 3. impose the common parameter structure and eliminate the excess parameters
%    the number of columns of the derivative matrices equals then the number
%    of FREE model parameters
DerivVecG = MIMO_ML_AddSelectColumns(DerivVecG, Sel, ModelVar);

% 4. reshape the derivative matrices to their original form
FreeParam.Theta = size(DerivVecG, 2);
DerivVecG = permute(reshape(DerivVecG, [nu*ny, F, FreeParam.Theta]), [1, 3, 2]);

% 5. calculation of the covariance matrices of the plant transfer function
CRbound.vecG = zeros(nu*ny, nu*ny, F);
CRbound.G = zeros(ny, nu, F);
for kk = 1:F
    XX = DerivVecG(:, :, kk) * SqrtCRtheta;
    CRbound.vecG(:, :, kk) = XX * XX';
    CRbound.G(:, :, kk) = reshape(diag(CRbound.vecG(:, :, kk)), [ny, nu]);
end % kk

