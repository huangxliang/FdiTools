function [CRbound, sqrtCRbound, TheCond] = MIMO_ML_CRtheta(PolyTrans, Deriv, Sel, ModelVar, data);
%
% function [CRbound, sqrtCRbound, TheCond] = MIMO_ML_CRtheta(PolyTrans, Deriv, Sel, ModelVar, data);
%
%
%   Output parameter
%
%		CRbound		=	Cramer-Rao bound of the estimated model parameters, the estimated plant model, and the estimated noise model
%							CRbound = struct('A', [], 'vecB',[], 'AvecB', [], 'all', [])
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
%                                   s-, sqrt(s)-domains:    a(OrderA) = 1
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
%							PolyTrans.A             =	denominator polynomial plant transfer function evaluated in x.Plant 
%                                                       size 1 x F 
%							PolyTrans.G             =	plant transfer matrix evaluated in x.Plant
%                                                       size ny x nu x F 
%							PolyTrans.Tg            =	plant transient term evaluated in x.Plant
%                                                       size ny x F 
%                           PolyTrans.sqrtCEinv     =   hermitian symmetric square root of the inverse of the covariance of the 
%                                                       output error (Cov(NY-G*NU)) 
%
%		Deriv		=	structure containing the derivative of vec(G), vec(H), vec(G'), vec(H') w.r.t.
%						all the plant model parameters a, b, c, d
%						    Deriv.vecGa             =	derivative vec(G) w.r.t. a;	size ny*nu x (na+1) x F 
%						    Deriv.vecGb             =	derivative vec(G) w.r.t. b;	size ny*nu x ny*nu*(nb+1) x F
%						    Deriv.vecGHa            =	derivative vec(G') w.r.t. a; size ny*nu x (na+1) x F
%						    Deriv.vecGHb            =	derivative vec(G') w.r.t. b; size ny*nu x ny*nu*(nb+1) x F
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
%															'':     plane not defined
%							ModelVar.Struct			=	model structure
%                                                           'EIV':	errors-in-variables (noisy input-output data)
%                                                           'OE':	generalised output error (known input, noisy output)
%							ModelVar.RecipPlant		=	1 if plant model is reciprocal: G(i,j) = G(j,i)
%							ModelVar.nu				=	number of inputs
%							ModelVar.ny				= 	number of outputs
%							ModelVar.na				=	order polynomial A
%							ModelVar.nb				= 	order matrix polynomial B
%
%		data		=	structure containing the non-parametric data required for the identification
%									data.U          =	For random signals one of the two following cases: 
%                                                           - input DFT spectrum 1 MIMO experiment:     nu x F 
%                                                           - power spectrum 1 MIMO experiment:         nu x nu x F 
%                                                       For periodic signals one of the two following cases: 
%                                                           - input DFT spectrum of 1 MIMO experiment:	nu x F 
%                                                           - input DFT specta of nu MIMO experiments:	nu x nu x F 
%									data.freq       =	vector of frequency values (Hz), size: F x 1
%									data.Ts         =	sampling time (s)
%									data.CY         =	(sample) noise covariance matrix of Y.
%                                                       For random signals:
%                                                           - covariance of 1 MIMO experiment:  ny x ny x F 
%                                                       For periodic signals one of the two following cases:  
%                                                           - 1 MIMO experiment:                ny x ny x F 
%                                                           - nu MIMO experiments:              ny x ny x nu x F 
%                                   data.CU         =   (sample) noise covariance matrix of U  
%                                                       For random signals:
%                                                           - covariance of 1 MIMO experiment:  nu x nu x F 
%                                                       For periodic signals one of the two following cases:  
%                                                           - 1 MIMO experiment:                nu x nu x F 
%                                                           - nu MIMO experiments:              nu x nu x nu x F 
%                                   data.CYU        =   (sample) noise covariance matrix of U 
%                                                       For random signals:
%                                                           - covariance of 1 MIMO experiment:  ny x nu x F 
%                                                       For periodic signals one of the two following cases:  
%                                                           - 1 MIMO experiment:                ny x nu x F 
%                                                           - nu MIMO experiments:              ny x nu x nu x F 
%                                   data.dof        =   degrees of freedom of the sample covariance estimates
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 24 October 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

na = ModelVar.na;                               % order polynomial A
nb = ModelVar.nb;                               % order matrix polynomial B
nig = ModelVar.nig;                             % order vector polynomial Ig
nu = ModelVar.nu;                               % number of inputs
ny = ModelVar.ny;                               % number of outputs
ntheta = (na+1) + (nb+1)*nu*ny + (nig+1)*ny;	% total number of model parameters
F = length(data.freq);                          % number of frequencies 

% check the number of MIMO experiments via the size 
% of the noise covariance matrix
sizeCov = size(data.CY);
if length(sizeCov) == 4
    NumberExp = size(data.CY, 3);               % nexp MIMO experiments
else
    NumberExp = 1;                              % 1 MIMO experiment
end % if

% intermediate data structure for calculating the CR-bound
dataee = struct('CY', [], 'CU', [], 'CYU', []);
dataee.CY = zeros(ny, ny, F);
switch ModelVar.Struct
    case 'EIV'
        dataee.CU = zeros(nu, nu, F);           % input covariance of 1 MIMO experiment
        dataee.CYU = zeros(ny, nu, F);          % output-input covariance of 1 MIMO experiment
    case 'OE'
        dataee.sqrtCYinv = zeros(ny, ny, F);
end % switch

CRbound = struct('A', [], 'AvecB', [], 'vecB', [], 'all', []);
% free parameters in each (matrix or vector) polynomial
FreeParam.A = sum(Sel.A);
FreeParam.B = sum(sum(sum(Sel.B)));
FreeParam.Theta = FreeParam.A + FreeParam.B;
CRbound.A = zeros(FreeParam.A, FreeParam.A);
CRbound.AvecB = zeros(FreeParam.A, FreeParam.B);
CRbound.vecB = zeros(FreeParam.B, FreeParam.B);
CRbound.all = zeros(FreeParam.Theta, FreeParam.Theta);

% input (power) spectrum
DimU = length(size(data.U));    % dimension matrix data.U
switch DimU
    case 2,
        LL = 1;                 % DFT spectrum nu x 1 reference signal is given, size nu x 1 x F
    case 3,
        if (NumberExp == 1) && (length(sizeCov) == 3)
            LL = nu;            % power spectrum nu x 1 reference signal is given, size nu x nu x F
        else % more than 1 MIMO experiment
            LL = 1;             % 1 MIMO experiment at the time should be handled
        end % if 1 MIMO experiment
end

% calculates the square root input power spectrum for random inputs only
% data.U = nu x nu x F matrix and data.CY = ny x ny x F matrix (1 experiment with a random signal) 
if (DimU == 3) && (length(sizeCov) == 3)
    % calculate a hermitian symmetric square root of data.U
    UR = zeros(nu, nu, F);
    for kk = 1:F
        [uR, sR, vR] = svd(data.U(:,:,kk),0);           % data.U is a hermitian symmetric positive definite matrix
        SqrtSR = uR * diag(diag(sR).^(0.5)) * uR';      % hermitian symmetric square root
        UR(:, :, kk) = SqrtSR;
    end % kk, frequencies
else % periodic signal
    UR = zeros(nu, F);
end % if input power spectrum is provided

% Jacobian matrix for all MIMO experiments
% If LL = nu then NumberExp = 1
TheJacob = zeros(LL*NumberExp*ny*F, ntheta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobian matrix w.r.t. all model parameters and all MIMO experiments %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ee = 1:NumberExp
    
    % data of MIMO experiment no. ee
    if length(sizeCov) == 4 % periodic excitations only
        UR(:,:) = data.U(:,ee,:);
        dataee.CY(:,:,:) = data.CY(:,:,ee,:);
        if strcmp(ModelVar.Struct, 'EIV')
            dataee.CU(:,:,:) = data.CU(:,:,ee,:);
            dataee.CYU(:,:,:) = data.CYU(:,:,ee,:);              
        end % if errors-in-variables
    else % 1 MIMO experiment; random or periodic excitations
        if LL == 1
            UR(:,:) = data.U;
        end % if not the input power spectrum is provided as input data
        dataee.CY(:,:,:) = data.CY;
         if strcmp(ModelVar.Struct, 'EIV')
            dataee.CU(:,:,:) = data.CU;
            dataee.CYU(:,:,:) = data.CYU;  
         end % if errors-in-variables
    end % if more than 1 MIMO experiment    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate a hermitian square root of the inverse %
    % of the covariance of the output error (NY-G*NU)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch ModelVar.Struct

        case 'EIV'
            PolyTrans = MIMO_ML_InvCovOutputError(dataee, PolyTrans);

        case 'OE'
            % relative upper limit under which the singular values are set to zero
            % in the pseudo-inverse of the square root of the weighting matrix 
            delta = 1e-6;
            % pseudo-inverse square root covariance matrix
            PolyTrans.sqrtCEinv = Sqrt_Inv(dataee.CY, delta);

    end % switch

    % inverse of a hermitian symmetric square root 
    % of the covariance of the output error (NY-G*NU) 
    sqrtCYinv = PolyTrans.sqrtCEinv;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Jacobian vec(W) w.r.t. all plant model parameters where          %
    % W = CY^(-1/2) * G  * U                                           %
    % size W:                                                          %
    %              ny x F for LL = 1                                   %
    %              ny x nu x F for LL = nu                             %
    % JacW = Jacobian vec(W) w.r.t. Theta, dimension ny x ntheta x F   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    JacW = zeros(LL*ny, ntheta, F);     

    for kk = 1:F

        % derivatives w.r.t. a
        Low = 1;
        Upp = (na+1);
        switch DimU
            case 2,
                MatW = kron(UR(:, kk).', sqrtCYinv(:,:,kk));
            case 3,
                if (NumberExp == 1) && (length(sizeCov) == 3) % random input, input power spectrum provided 
                    MatW = kron(UR(:,:,kk).', sqrtCYinv(:,:,kk)); 
                else % periodic input, nu MIMO experiments
                    MatW = kron(UR(:,kk).', sqrtCYinv(:,:,kk));                    
                end % if 1 MIMO experiment
        end
        JacW(:, Low:Upp, kk) = MatW * Deriv.vecGa(:,:,kk);

        % derivatives w.r.t. b
        Low = Upp + 1;
        Upp = Low + (nb+1)*ny*nu - 1;
        JacW(:, Low:Upp, kk) = MatW * Deriv.vecGb(:,:,kk);

    end % frequencies kk

    % DC and Nyquist count for 1/2 in the sums => as 1/sqrt(2) in the Jacobian matrices
    if data.DC == 1
        JacW(:,:,1) = JacW(:,:,1)/sqrt(2);
    end % if DC
    if data.Nyquist == 1
        JacW(:,:,end) = JacW(:,:,end)/sqrt(2);
    end % if Nyquist


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Put the Jacobian in the following size:                           %
    %       JacW: LL*ny*F x ntheta                                         %
    % 2. impose common parameter structure and eliminate excess parameters %
    % 3. put real and imaginary parts on top of each other                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % put all frequency contributions on top of each other
    JacW = reshape(permute(JacW, [1,3,2]), LL*ny*F, ntheta);
    
    % put experiment no. ee in the Jacobian matrix and prediction error vector
    SelectRows = [(ee-1)*LL*ny*F+1:ee*LL*ny*F];
    TheJacob(SelectRows, :) = JacW;
    
end % ee, MIMO experiments

% impose common parameter structure and eliminate excess parameters
TheJacob = MIMO_ML_AddSelectColumns(TheJacob, Sel, ModelVar);

% put real and imaginary parts on top of each other
% factor sqrt(2) needed because the noise cov. of the complex values is used 
TheJacob = sqrt(2)*[real(TheJacob); imag(TheJacob)]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CR-bound free parameters Theta %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[uj, sj, vj] = svd(TheJacob, 0);

% final CR-bound estimated model parameters
vj = vj * diag(1./diag(sj));
sqrtCRbound = vj;
CRbound.all = vj * vj.';
TheCond = sj(1,1)/sj(end,end);

% scaling CRbound for the variability of the estimated noise covariances
if isinf(data.dof)
    ScaleCRML = 1;                              % exactly known covariances
else
    dof = data.dof;
    ScaleCRML = (dof)^2/(dof+1-ny)/(dof-ny-1);  % estimated covariances
end % if
CRbound.all = ScaleCRML * CRbound.all;
sqrtCRbound = sqrt(ScaleCRML) * sqrtCRbound;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CR-bound free parameters (matrix or vector) polynomials %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% free parameters A-polynomial
Low = 1;
Upp = Low + FreeParam.A - 1;
CRbound.A = CRbound.all(Low:Upp, Low:Upp);       % CR-bound
LowA = Low; UppA = Upp;

% free parameters B-polynomial matrix
Low = Upp + 1;
Upp = Low + FreeParam.B - 1;
CRbound.vecB = CRbound.all(Low:Upp, Low:Upp);         % CR-bound

% covariance between free A- and B-coefficients
CRbound.AvecB = CRbound.all(LowA:UppA, Low:Upp);      % CR-bound

