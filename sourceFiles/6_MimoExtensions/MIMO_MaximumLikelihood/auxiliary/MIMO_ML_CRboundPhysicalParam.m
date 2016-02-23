function CRbound = MIMO_ML_CRboundPhysicalParam(CRbound, Seln, wscale, ModelVar);
%
%       CRbound = MIMO_ML_CRboundPhysicalParam(CRbound, Seln, wscale, ModelVar);
%
%
%	Output parameters
%
%		CRbound				=	Cramer-Rao bound of the estimated physical model parameters, and the estimated plant model 
%								CRbound = struct('A', [], 'AvecB', [], 'vecB', [], 'res', [], 'poles', [])
%                               See Input parameters for the definition of the structure 
%
%
%	Input parameters
%
%		CRbound				=	Cramer-Rao bound of the estimated normalised model parameters, the estimated plant model, and the estimated noise model
%								CRbound = struct('A', [], 'AcecB', [], 'vecB', [], 'res', [], 'poles', [])
%									CRbound.A               =   FreeParam.A x FreeParam.A
%                                                                   CRbound.A(i,j)      = covariance between coefficients a(i-1) 
%                                                                                         and a(j-1) 
%							        CRbound.AvecB           =   FreeParam.A x FreeParam.B
%                                                                   CRbound.AvecB(i,j)	= covariance between free coefficients 
%                                                                                         a(i-1) and vecB(j) 
%									CRbound.vecB            =   FreeParam.B x FreeParam.B
%                                                                   CRbound.vecB(i,j)   = covariance between vecB(i) and vecB(j) 
%                                                                                         where vecB = permute(B, [3, 1, 2]); 
%                                                                                         vecB = vecB(:)
%                                   CRbound.res             =   struct{'all', 'sv', 'lsv', 'rsv'}
%                                                                   CRbound.res.all     =   covariance matrix of (vec(Res))re; where ()re puts the real
%                                                                                           and imaginary parts of the matrix on top of each other; vec()
%                                                                                           stacks the columns of the matrix on top of each other; and Res
%                                                                                           is a residue matrix of the ny x nu transfer function matrix G = B/A;
%                                                                                           size: (2*ny*nu) x (2*ny*nu) x na
%                                                                   CRbound.res.sv      =   variance singular values of the residues;
%                                                                                           size: min(ny, nu) x na
%                                                                   CRbound.res.lsv     =   covariance left singular vectors of the residues [real(ur); imag(ur)];
%                                                                                           size: 2*ny x 2*ny x min(ny, nu) x na
%                                                                   CRbound.res.rsv     =   covariance right singular vectors of the residues [real(vr); imag(vr)];
%                                                                                           size: 2*nu x 2*nu x min(ny, nu) x na
%                                   CRbound.poles           =   struct{'root', 'all', 'damp', 'freq', 'time'}
%                                                                   CRbound.poles.root	=   cov((root)re); where ()re stacks the real and imaginary part of the root
%                                                                                           on top of each other; ; size 2 x 2 x number of roots
%                                                                   CRbound.poles.all	=   Cov(roots.all) = cov((roots)re); where the real and imaginary parts of
%                                                                                           the vector roots are stacked on top of each other;
%                                                                                           size 2*(number of roots) x 2*(number of roots)
%                                                                   CRbound.poles.damp	=   variance damping complex roots; entry is NaN for real roots
%                                                                   CRbound.poles.freq	=   variance frequency complex roots; entry is NaN for real roots
%                                                                   CRbound.poles.time	=   variance time constant real roots; entry is NaN for complex roots
%
%		Seln				=	struct('A', [], 'B', [])
%									Seln.A       =   1 x (OrderA+1)
%                                                       Seln.A(r) = 1 if coeff. a(r-1) is unknown
%                                                       Seln.A(r) = 0 if coeff. a(r-1) = 0
%									Seln.B       =   ny x nu x (OrderB+1)
%                                                       Seln.B(i,j,r) = 1 if coeff. b(i,j,r-1) is unknown
%                                                       Seln.B(i,j,r) = 0 if coeff. b(i,j,r-1) = 0
%
%		wscale				=	angular the frequency scaling
%
%		ModelVar			=	contains the information about the model to be identified
%								struct('Transient', [], 'PlantPlane', [], 'Struct', [], 'RecipPlant', [])
%									ModelVar.Transient		=	1 then the initial conditions of the plant are estimated
%									ModelVar.PlantPlane		=	plane of the plant model
%																	's':	continuous-time;
%																	'w':	sqrt(s)-domain
%																	'z':	discrete-time;
%																	'':		plane not defined
%									ModelVar.Struct			=	model structure
%                                                                   'EIV':  errors-in-variables (noisy input-output data)
%                                                                   'OE':	generalised output error (known input, noisy output)
%									ModelVar.RecipPlant		=	1 if plant model is reciprocal: G(i,j) = G(j,i)
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 6 October 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plant model parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(ModelVar.PlantPlane, 'z')
    
    na = size(Seln.A, 2) - 1;
    nb = size(Seln.B, 3) - 1;
    nu = size(Seln.B, 2);
    ny = size(Seln.B, 1);
    nn = max(na, nb) + 1;
    Scale = ones(nn, 1);
    for ii = 1:nn
        Scale(ii) = wscale^(ii-1);
    end % for ii
    Scale = Scale * Scale.';
    
    % scaling CR-bound coefficients matrix polynomial B
    nb_all = sum(sum(sum(Seln.B)));
    ScaleMat  = zeros(nb_all, nb_all);
    OffsetRow = 0;
    for ll = 1:nu
        for kk = 1:ny
            nrow = sum(Seln.B(kk,ll,:));
            OffsetColumn = 0;
            for jj = 1:nu
                for ii = 1:ny
                    ncolumn = sum(Seln.B(ii,jj,:));
                    ScaleMat(OffsetRow+1:OffsetRow+nrow, OffsetColumn+1:OffsetColumn+ncolumn) = ...
                             Scale(Seln.B(kk,ll,:) == 1, Seln.B(ii,jj,:) == 1);
                    OffsetColumn = OffsetColumn + ncolumn;
                end % ii row index
            end % jj column index
            OffsetRow = OffsetRow + nrow;
        end % kk row index
    end %ll column index
    CRbound.vecB = CRbound.vecB ./ ScaleMat;
    
    % scaling CR-bound coefficients A-polynomial
    CRbound.A = CRbound.A ./ Scale(Seln.A==1, Seln.A==1);
    
    % scaling CR-bound covariance A- and B-coefficients
    na_all = sum(Seln.A);
    ScaleMat = zeros(na_all, nb_all);
    OffsetColumn = 0;
    for jj = 1:nu
        for ii = 1:ny
            ncolumn = sum(Seln.B(ii,jj,:));
            ScaleMat(:, OffsetColumn+1:OffsetColumn+ncolumn) = Scale(Seln.A == 1, Seln.B(ii,jj,:) == 1);
            OffsetColumn = OffsetColumn + ncolumn;
        end % ii row index
    end % jj column index
    CRbound.AvecB = CRbound.AvecB ./ ScaleMat; 
     
    switch ModelVar.PlantPlane 
        case 's'
            wscale2 = wscale;
        case 'w'
            % poles are squared in sqrt(s)-domain before calculating
            % the resonance frequencies and damping ratios
            wscale2 = wscale^2;   
    end % if
     
    % scaling CR-bound of the poles, frequency, and time
    % constants, but not the damping because it has no units
    CRbound.poles.all = CRbound.poles.all * (wscale^2);
    CRbound.poles.root = CRbound.poles.root * (wscale^2);
    CRbound.poles.freq = CRbound.poles.freq * (wscale2^2);
    CRbound.poles.time = CRbound.poles.time / (wscale2^2);

    % scaling CR-bound of the residue matrices and the singular values, 
    % but not the left and right singular vectors because they have no units
    CRbound.res.all = CRbound.res.all * (wscale^2);
    CRbound.res.sv = CRbound.res.sv * (wscale^2);
   
end % if not z-domain

