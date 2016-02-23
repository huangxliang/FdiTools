function [CRbound, sqrtCRbound, TheCond] = CRmodelparameters(PolyTrans, Deriv, Sel, ModelVar, data);
%
% function [CRbound, sqrtCRbound, TheCond] = CRmodelparameters(PolyTrans, Deriv, Sel, ModelVar, data);
%
%
%   Output parameter
%
%		CRbound		=	Cramer-Rao bound of the estimated model parameters, the estimated plant model, and the estimated noise model
%						structure with fields 'A', 'vecB', 'vecC', 'D', 'Theta'
%							CRbound = struct('A', [], 'vecB',[], 'vecC', [], 'D', [], 'Theta', [])
%							CRbound.A = FreeParam.A x FreeParam.A
%								CRbound.A(i,j) = covariance between free coefficients a(i-1) and a(j-1) 
%							CRbound.AvecB = FreeParam.A x FreeParam.B
%								CRbound.AvecB(i,j) = covariance between free coefficients a(i-1) and vecB(j) 
%							CRbound.vecB = FreeParam.B x FreeParam.B
%								CRbound.vecB(i,j) = covariance between free parameters vecB(i) and vecB(j)
%							CRbound.vecC = FreeParam.C x FreeParam.C
%								CRbound.vecC(i,j) = covariance between free parameters vecC(i) and vecC(j) 
%							CRbound.D = FreeParam.D x FreeParam.D
%								CRbound.D(i,j) = covariance between free coefficients d(i-1) and d(j-1) 
%							CRbound.DvecC = FreeParam.D x FreeParam.C
%								CRbound.DvecC(i,j) = covariance between free coefficients d(i-1) and vecC(j) 
%							CRbound.all = dim(Theta) x dim(Theta), where dim(Theta) = FreeParam.A + FreeParam.B + FreeParam.C + FreeParam.D
%								CRbound.all(i,j) = covariance between free parameters Theta(i) and Theta(j)
%							Notes:
%                               - in s-, sqrt(s) domains the CR-bound of the normalised parameters is calculated
%                               - the (normalised) model parameters satisfy the following constraints
%                                   z-domain:               a(0) = d(0) = 1
%                                   s-, sqrt(s)-domains:    a(OrderA) = d(OrderD) = 1
%                               - the C-coefficients are normalised such that CovE = eye(ny)
%
%       sqrtCRbound =   square root of Cramer-Rao lower bound on all the model parameters to guarantee a numerical stable calculation of
%                           dim(Theta) x dim(Theta); and CRbound.all = sqrtCRbound * sqrtCRbound.'
%
%       TheCond     =   condition number Jacobian matrix used to calculate the Cramer-Rao bound (cond(CRbound) = TheCond^2)
%
%
%   Input parameters
%
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
%		data		=	structure containing the non-parametric data required for the identification
%							data.Y		    =	DFT spectrum ny x 1 output signal, dimensions ny x number of frequencies
%							data.U		    =	DFT spectrum nu x 1 input signal, dimensions: nu x number of frequencies
%                           data.R          =   either
%                                                   - DFT spectrum nu x 1 reference signal, dimensions: nu x number of frequencies
%                                                   - power spectrum nu x 1 reference signal, dimensions: nu x nu x number of frequencies
%							data.freq	    =	vector of frequency values (Hz), dimension: number of frequencies x 1
%							data.Ts		    =	sampling time (s)
%							data.Gc		    =	controller transfer function, dimension nu x ny x number of frequencies
%                           data.DC         =   1 if DC is in the frequency set
%                           data.Nyquist    =   1 if Nyquist is in the frequency set
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, April 2005 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 18 October 2011
%

na = ModelVar.na;
nb = ModelVar.nb;
nc = ModelVar.nc;
nd = ModelVar.nd;
nig = ModelVar.nig;
nih = ModelVar.nih;
nu = ModelVar.nu;
ny = ModelVar.ny;
ntheta = (na+1) + (nb+1)*nu*ny + (nig+1)*ny + (nc+1)*ny^2 + (nd+1) + (nih+1)*ny; % total number of model parameters
F = length(data.freq);

CRbound = struct('A', [], 'vecB',[], 'vecC', [], 'D', [], 'all', []);
% free parameters in each (matrix or vector) polynomial
FreeParam.A = sum(Sel.A);
FreeParam.B = sum(sum(sum(Sel.B)));
FreeParam.C = sum(sum(sum(Sel.C)));
FreeParam.D = sum(Sel.D);
FreeParam.Theta = FreeParam.A + FreeParam.B + FreeParam.C + FreeParam.D;
CRbound.A = zeros(FreeParam.A, FreeParam.A);
CRbound.vecB = zeros(FreeParam.B, FreeParam.B);
CRbound.vecC = zeros(FreeParam.C, FreeParam.C);
CRbound.D = zeros(FreeParam.D, FreeParam.D);
CRbound.all = zeros(FreeParam.Theta, FreeParam.Theta);

% intermediate variables used to calculate the derivatives
Gc_S = zeros(nu, ny, F);
UR = zeros(size(data.R));
DimR = length(size(UR));        % dimension matrix data.R
switch DimR
    case 2,
        LL = 1;                 % DFT spectrum ny x 1 reference signal is given, size ny x 1 x F
    case 3,
        LL = nu;                % power spectrum ny x 1 reference signal is given, size ny x ny x F
end
for kk = 1:F
    Gc_S(:, :, kk) = data.Gc(:, :, kk) * inv(eye(ny) + PolyTrans.G(:, :, kk) * data.Gc(:, :, kk)) * PolyTrans.H(:, :, kk);
    switch DimR
        case 2,
            UR(:, kk) = inv(eye(nu) + data.Gc(:, :, kk) * PolyTrans.G(:, :, kk)) * data.R(:, kk);
        case 3,
            [uR, sR, vR] = svd(data.R(:,:,kk),0);           % data.R is a hermitian symmetric positive definite matrix
            SqrtSR = uR * diag(diag(sR).^(0.5)) * uR';      % hermitian symmetric square root
            UR(:, :, kk) = inv(eye(nu) + data.Gc(:, :, kk) * PolyTrans.G(:, :, kk)) * SqrtSR;
    end  
end % kk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobian vec(W) and vec(V) w.r.t. all model parameters where     %
% W = H0^(-1) * Delta_G * [eye(nu) + Gc * G0]^(-1) * R             %
% dimension W:                                                     %
%              ny x F for LL = 1                                   %
%              ny x nu x F for LL = nu                             %
% JacW = Jacobian vec(W) w.r.t. Theta, dimension ny x ntheta x F   %
% V =  S0^(-1) * S * S' * S0'^(-1)                                 %
% dimension V: ny x ny x F                                         %
% JacV = Jacobian vec(V) w.r.t. Theta, dimension ny^2 x ntheta x F %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JacW = zeros(LL*ny, ntheta, F);     
JacV = zeros(ny^2, ntheta, F);

for kk = 1:F
	
	% derivatives w.r.t. a
	Low = 1;
	Upp = (na+1);
    switch DimR
        case 2,
	        MatW = kron(UR(:, kk).', PolyTrans.Hinv(: ,: , kk));
        case 3,
	        MatW = kron(UR(:, :, kk).', PolyTrans.Hinv(: ,: , kk));            
    end
	JacW(:, Low:Upp, kk) = MatW * Deriv.vecGa(: ,: , kk);
    MatV1 = kron(Gc_S(:, :, kk).', PolyTrans.Hinv(:, :, kk));
    MatV2 = kron(conj(PolyTrans.Hinv(:, :, kk)), Gc_S(:, :, kk)');
    JacV(:, Low:Upp, kk) = - MatV1 * Deriv.vecGa(:, :, kk) - MatV2 * Deriv.vecGHa(:, :, kk);

	% derivatives w.r.t. b
	Low = Upp + 1;
	Upp = Low + (nb+1)*ny*nu - 1;
	JacW(:, Low:Upp, kk) = MatW * Deriv.vecGb(: ,: , kk);
    JacV(:, Low:Upp, kk) = - MatV1 * Deriv.vecGb(:, :, kk) - MatV2 * Deriv.vecGHb(:, :, kk);
	
	% derivatives w.r.t. ig are zero
	Low = Upp + 1;
	Upp = Low + (nig+1)*ny - 1;
	
	% derivatives w.r.t. c
	Low = Upp + 1;
	Upp = Low + (nc+1)*ny^2 - 1;
    MatV1 = kron(eye(ny), PolyTrans.Hinv(:, :, kk));
    MatV2 = kron(conj(PolyTrans.Hinv(:, :, kk)), eye(ny));
	JacV(:, Low:Upp, kk) = MatV1 * Deriv.vecHc(:, :, kk) + MatV2 * Deriv.vecHHc(:, :, kk);

	% derivatives w.r.t. d
	Low = Upp + 1;
	Upp = Low + (nd+1) - 1;
	JacV(:, Low:Upp, kk) = MatV1 * Deriv.vecHd(:, :, kk) + MatV2 * Deriv.vecHHd(:, :, kk);
	
% 	% derivatives w.r.t. ih are zero
% 	Low = Upp + 1;
% 	Upp = Low + (nih+1)*ny - 1;
	
end % frequencies kk

% DC and Nyquist count for 1/2 in the sums => as 1/sqrt(2) in the Jacobian matrices
F1 = F;
if data.DC == 1
    F1 = F1-1/2;
    JacW(:,:,1) = JacW(:,:,1)/sqrt(2);
    JacV(:,:,1) = JacV(:,:,1)/sqrt(2);
end % if DC
if data.Nyquist == 1
    F1 = F1-1/2;
    JacW(:,:,end) = JacW(:,:,end)/sqrt(2);
    JacV(:,:,end) = JacV(:,:,end)/sqrt(2);
end % if Nyquist


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Put Jacobians in the following size:                              %
%       JacW: LL*ny*F x ntheta                                         %
%       JacV: ny^2*F x ntheta                                          %
% 2. impose common parameter structure and eliminate excess parameters %
% 3. put real and imaginary parts on top of each other                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put all frequency contributions on top of each other
JacW = reshape(permute(JacW, [1,3,2]), LL*ny*F, ntheta);
JacV = reshape(permute(JacV, [1,3,2]), ny^2*F, ntheta);

% impose common parameter structure and eliminate excess parameters
JacW = Add_SelectColumns(JacW, Sel, ModelVar);
JacV = Add_SelectColumns(JacV, Sel, ModelVar);

% put real and imaginary parts on top of each other
JacW = [real(JacW); imag(JacW)];                % size 2*ny*F x ntheta
JacV = [real(JacV); imag(JacV)];                % size 2*ny^2*F x ntheta


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CR-bound free parameters Theta %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full Jacobian matrix
Jacob = [sqrt(2)*JacW; JacV];
[uj, sj, vj] = svd(Jacob, 0);

% final CR-bound estimated model parameters
vj = vj * diag(1./diag(sj));
sqrtCRbound = vj;
CRbound.all = vj * vj.';
TheCond = sj(1,1)/sj(end,end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CR-bound free parameters (matrix or vector) polynomials %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% free parameters A-polynomial
Low = 1;
Upp = Low + FreeParam.A - 1;
CRbound.A = CRbound.all(Low:Upp, Low:Upp);              % CR-bound
LowA = Low; UppA = Upp;

% free parameters B-polynomial matrix
Low = Upp + 1;
Upp = Low + FreeParam.B - 1;
CRbound.vecB = CRbound.all(Low:Upp, Low:Upp);           % CR-bound

% covariance between free A- and B-coefficients
CRbound.AvecB = CRbound.all(LowA:UppA, Low:Upp);        % CR-bound

% free parameters C-polynomial matrix
Low = Upp + 1;
Upp = Low + FreeParam.C - 1;
CRbound.vecC = CRbound.all(Low:Upp, Low:Upp);           % CR-bound
LowC = Low; UppC = Upp;

if ~strcmp(ModelVar.Struct, 'ARMAX')
    
    % free parameters D-polynomial
    Low = Upp + 1;
    Upp = Low + FreeParam.D - 1;
    CRbound.D = CRbound.all(Low:Upp, Low:Upp);          % CR-bound

    % covariance between free D- and C-coefficients
    CRbound.DvecC = CRbound.all(Low:Upp, LowC:UppC);	% CR-bound 
    
else % ARMAX model
    
    % noise denominator = plant denominator
    CRbound.D = CRbound.A; 
    CRbound.DvecC = CRbound.all(LowA:UppA, LowC:UppC);

end % if not ARMAX model
