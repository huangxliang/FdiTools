function CRbound = CRboundPhysicalParam(CRbound, Seln, wscale, ModelVar);
%
%       CRbound = CRboundPhysicalParam(CRbound, Seln, wscale, ModelVar);
%
%
%	Output parameters
%
%		CRbound				=	Cramer-Rao bound of the estimated physical model parameters, the estimated plant model, and the estimated noise model
%								CRbound = struct('A', [], 'vecB',[], 'vecC', [], 'D', [], 'plant', [], 'noise', []); 
%                               definition struct: see input parameters 
%
%
%	Input parameters
%
%		CRbound				=	Cramer-Rao bound of the estimated normalised model parameters, the estimated plant model, and the estimated noise model
%								CRbound = struct('A', [], 'vecB',[], 'vecC', [], 'D', [], 'plant', [], 'noise', [])
%									CRbound.A               =   FreeParam.A x FreeParam.A
%                                                                   CRbound.A(i,j) = covariance between coefficients a(i-1) and a(j-1) 
%							        CRbound.AvecB           =   FreeParam.A x FreeParam.B
%                                                                   CRbound.AvecB(i,j) = covariance between free coefficients a(i-1) and vecB(j) 
%									CRbound.vecB            =   FreeParam.B x FreeParam.B
%                                                                   CRbound.vecB(i,j) = covariance between vecB(i) and vecB(j) where vecB = permute(B, [3, 1, 2]); vecB = vecB(:)
%									CRbound.vecC            =   FreeParam.C x FreeParam.C
%                                                                   CRbound.vecC(i,j) = covariance between vecC(i) and vecC(j) where vecC = permute(C, [3, 1, 2]); vecC = vecC(:) 
%									CRbound.D               =   FreeParam.D x FreeParam.D
%                                                                   CRbound.D(i,j) = covariance between coefficients d(i-1) and d(j-1) 
%							        CRbound.DvecC           =   FreeParam.D x FreeParam.C
%                                                                   CRbound.DvecC(i,j) = covariance between free coefficients d(i-1) and vecC(j) 
%                                   CRbound.plant           =   struct('res', [], 'poles', []) 
%                                                                   covariance of the plant residue matrices and the plant poles  
%                                   CRbound.noise.res       =   struct('all', 'sv', 'lsv', 'rsv') 
%                                                                   CRbound.noise.res.all	=   covariance matrix of (vec(Res))re; where ()re puts the real
%                                                                                               and imaginary parts of the matrix on top of each other; vec()
%                                                                                               stacks the columns of the matrix on top of each other; and Res
%                                                                                               is a residue matrix of the ny x ny noise transfer function matrix H = C/D;
%                                                                                                   size: (2*ny^2) x (2*ny^2) x na
%                                                                   CRbound.noise.res.sv	=   variance singular values of the residues;
%                                                                                                   size: ny x na
%                                                                   CRbound.noise.res.lsv	=   covariance left singular vectors of the residues [real(ur); imag(ur)];
%                                                                                                   size: 2*ny x 2*ny x ny x na
%                                                                   CRbound.noise.res.rsv   =   covariance right singular vectors of the residues [real(vr); imag(vr)];
%                                                                                                   size: 2*ny x 2*ny x ny x na
%                                   CRbound.noise.poles     =   struct{'root', 'all', 'damp', 'freq', 'time'}
%                                                                   CRbound.noise.poles.root	=   cov((root)re); where ()re stacks the real and imaginary part of the root
%                                                                                                   on top of each other; ; size 2 x 2 x number of roots
%                                                                   CRbound.noise.poles.all     =   Cov(roots.all) = cov((roots)re); where the real and imaginary parts of
%                                                                                                   the vector roots are stacked on top of each other;
%                                                                                                       size 2*(number of roots) x 2*(number of roots)
%                                                                   CRbound.noise.poles.damp	=   variance damping complex roots; entry is NaN for real roots
%                                                                   CRbound.noise.poles.freq	=   variance frequency complex roots; entry is NaN for real roots
%                                                                   CRbound.noise.poles.time	=   variance time constant real roots; entry is NaN for complex roots
%                                   CRbound.plant.res       =   struct('all', 'sv', 'lsv', 'rsv') 
%                                                                   CRbound.plant.res.all	=   covariance matrix of (vec(Res))re; where ()re puts the real
%                                                                                               and imaginary parts of the matrix on top of each other; vec()
%                                                                                               stacks the columns of the matrix on top of each other; and Res
%                                                                                               is a residue matrix of the ny x nu plant transfer function matrix G = B/A;
%                                                                                                   size: (2*ny*nu) x (2*ny*nu) x na
%                                                                   CRbound.plant.res.sv	=   variance singular values of the residues;
%                                                                                                   size: min(ny, nu) x na
%                                                                   CRbound.plant.res.lsv	=   covariance left singular vectors of the residues [real(ur); imag(ur)];
%                                                                                                   size: 2*ny x 2*ny x min(ny, nu) x na
%                                                                   CRbound.plant.res.rsv   =   covariance right singular vectors of the residues [real(vr); imag(vr)];
%                                                                                                   size: 2*nu x 2*nu x min(ny, nu) x na
%                                   CRbound.plant.poles     =   struct{'root', 'all', 'damp', 'freq', 'time'}
%                                                                   CRbound.plant.poles.root	=   cov((root)re); where ()re stacks the real and imaginary part of the root
%                                                                                                   on top of each other; ; size 2 x 2 x number of roots
%                                                                   CRbound.plant.poles.all     =   Cov(roots.all) = cov((roots)re); where the real and imaginary parts of
%                                                                                                   the vector roots are stacked on top of each other;
%                                                                                                       size 2*(number of roots) x 2*(number of roots)
%                                                                   CRbound.plant.poles.damp	=   variance damping complex roots; entry is NaN for real roots
%                                                                   CRbound.plant.poles.freq	=   variance frequency complex roots; entry is NaN for real roots
%                                                                   CRbound.plant.poles.time	=   variance time constant real roots; entry is NaN for complex roots
%
%		Seln				=	struct('A', [], 'B', [], 'C', [], 'D', [])
%									Seln.A       =   1 x (OrderA+1)
%                                                       Seln.A(r) = 1 if coeff. a(r-1) is unknown
%                                                       Seln.A(r) = 0 if coeff. a(r-1) = 0
%									Seln.B       =   ny x nu x (OrderB+1)
%                                                       Seln.B(i,j,r) = 1 if coeff. b(i,j,r-1) is unknown
%                                                       Seln.B(i,j,r) = 0 if coeff. b(i,j,r-1) = 0
%									Seln.C       =   ny x ny x (OrderC+1)
%                                                       Seln.C(i,j,r) = 1 if coeff. c(i,j,r-1) is unknown
%                                                       Seln.C(i,j,r) = 0 if coeff. c(i,j,r-1) = 0
%									Seln.D       =   1 x (OrderD+1)
%                                                       Seln.D(r) = 1 if coeff. d(i,j,r-1) is unknown
%                                                       Seln.D(r) = 0 if coeff. d(i,j,r-1) = 0
%
%		wscale				=	structure containing the frequency scaling
%									wscale = struct('Plant', [], 'Noise', [])
%									wscale.Plant	=	angular frequency scaling plant model
%									wscale.Noise	=	angular frequency scaling noise model
%
%		ModelVar			=	contains the information about the model to be identified
%								struct('Transient', [], 'PlantPlane', [], 'NoisePlane', [], 'Struct', [], 'RecipPlant',[], 'RecipNoise')
%									ModelVar.Transient		=	1 then the initial conditions of the plant and/or noise are estimated
%									ModelVar.PlantPlane		=	plane of the plant model
%																	's':	continuous-time;
%																	'w':	sqrt(s)-domain
%																	'z':	discrete-time;
%																	'':		plane not defined
%									ModelVar.NoisePlane		=	plane of the plant model
%																	's':	continuous-time;
%																	'w':	sqrt(s)-domain
%																	'z':	discrete-time;
%																	'':		plane not defined
%									ModelVar.Struct			=	model structure
%																	'BJ':		Box-Jenkins
%																	'OE':		output error (plant model only)
%																	'ARMA':		autoregressive moving average (noise model only)
%																	'ARMAX':	autoregressive moving average with exogenous input
%									ModelVar.RecipPlant		=	1 if plant model is reciprocal: G(i,j) = G(j,i)
%									ModelVar.RecipNoise		=	1 if noise model is reciprocal: H(i,j) = H(j,i)
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, April 2008 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 18 October 2011
%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plant model parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(ModelVar.PlantPlane, 's') || strcmp(ModelVar.PlantPlane, 'w')
    
    na = size(Seln.A, 2) - 1;
    nb = size(Seln.B, 3) - 1;
    nu = size(Seln.B, 2);
    ny = size(Seln.B, 1);
    nn = max(na, nb) + 1;
    Scale = ones(nn, 1);
    for ii = 1:nn
        Scale(ii) = wscale.Plant^(ii-1);
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
    
    wscale1 = wscale.Plant;
    switch ModelVar.PlantPlane 
        case 's'
            wscale2 = wscale.Plant;
        case 'w'
            % poles are squared in sqrt(s)-domain before calculating
            % the resonance frequencies and damping ratios
            wscale2 = wscale.Plant^2;   
    end % if
     
    % scaling CR-bound of the plant poles, frequency, and time
    % constants, but not the damping because it has no units
    CRbound.plant.poles.all = CRbound.plant.poles.all * (wscale1^2);
    CRbound.plant.poles.root = CRbound.plant.poles.root * (wscale1^2);
    CRbound.plant.poles.freq = CRbound.plant.poles.freq * (wscale2^2);
    CRbound.plant.poles.time = CRbound.plant.poles.time / (wscale2^2);

    % scaling CR-bound of the plant residue matrices and the singular values, 
    % but not the left and right singular vectors because they have no units
    CRbound.plant.res.all = CRbound.plant.res.all * (wscale1^2);
    CRbound.plant.res.sv = CRbound.plant.res.sv * (wscale1^2);
    
end % if s- or sqrt(s)-domain plant model


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noise model parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(ModelVar.NoisePlane, 's') || strcmp(ModelVar.NoisePlane, 'w')
    
    nd = size(Seln.D, 2) - 1;
    nc = size(Seln.C, 3) - 1;
    ny = size(Seln.C, 1);
    nn = max(nc, nd) + 1;
    Scale = ones(nn, 1);
    for ii = 1:nn
        Scale(ii) = wscale.Noise^(ii-1);
    end % for ii
    Scale = Scale * Scale.';
    
    % scaling CR-bound coefficients matrix polynomial C
    nc_all = sum(sum(sum(Seln.C)));
    ScaleMat  = zeros(nc_all, nc_all);
    OffsetRow = 0;
    for ll = 1:ny
        for kk = 1:ny
            nrow = sum(Seln.C(kk,ll,:));
            OffsetColumn = 0;
            for jj = 1:ny
                for ii = 1:ny
                    ncolumn = sum(Seln.C(ii,jj,:));
                    ScaleMat(OffsetRow+1:OffsetRow+nrow, OffsetColumn+1:OffsetColumn+ncolumn) = ...
                             Scale(Seln.C(kk,ll,:) == 1, Seln.C(ii,jj,:) == 1);
                    OffsetColumn = OffsetColumn + ncolumn;
                end % ii row index
            end % jj column index
            OffsetRow = OffsetRow + nrow;
        end % kk row index
    end %ll column index
    CRbound.vecC = CRbound.vecC ./ ScaleMat;
    
    % scaling CR-bound coefficients D-polynomial
    CRbound.D = CRbound.D ./ Scale(Seln.D==1, Seln.D==1);
    
    % scaling CR-bound covariance D- and C-coefficients
    nd_all = sum(Seln.D);
    ScaleMat = zeros(nd_all, nc_all);
    OffsetColumn = 0;
    for jj = 1:ny
        for ii = 1:ny
            ncolumn = sum(Seln.C(ii,jj,:));
            ScaleMat(:, OffsetColumn+1:OffsetColumn+ncolumn) = Scale(Seln.D == 1, Seln.C(ii,jj,:) == 1);
            OffsetColumn = OffsetColumn + ncolumn;
        end % ii row index
    end % jj column index
    CRbound.DvecC = CRbound.DvecC ./ ScaleMat; 
         
    wscale1 = wscale.Noise;
    switch ModelVar.NoisePlane 
        case 's'
            wscale2 = wscale.Noise;
        case 'w'
            % poles are squared in sqrt(s)-domain before calculating
            % the resonance frequencies and damping ratios
            wscale2 = wscale.Noise^2;   
    end % if
     
    % scaling CR-bound of the noise poles, frequency, and time
    % constants, but not the damping because it has no units
    CRbound.noise.poles.all = CRbound.noise.poles.all * (wscale1^2);
    CRbound.noise.poles.root = CRbound.noise.poles.root * (wscale1^2);
    CRbound.noise.poles.freq = CRbound.noise.poles.freq * (wscale2^2);
    CRbound.noise.poles.time = CRbound.noise.poles.time / (wscale2^2);

    % scaling CR-bound of the noise residue matrices and the singular values, 
    % but not the left and right singular vectors because they have no units
    CRbound.noise.res.all = CRbound.noise.res.all * (wscale1^2);
    CRbound.noise.res.sv = CRbound.noise.res.sv * (wscale1^2);
    
end % if  s- or sqrt(s)-domain noise model
