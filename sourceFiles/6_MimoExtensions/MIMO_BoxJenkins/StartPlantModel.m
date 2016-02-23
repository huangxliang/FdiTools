function ThetaPlant = StartPlantModel(data, Sel, ModelVar, IterVar, FigNum)
%
%   ThetaPlant = StartPlantModel(data, Sel, ModelVar, IterVar, FigNum)
%
%       Generation of staring values for the plant model parameters G = B/A
%       and Tg = Ig/A using a weighted nonlinear least squares estimator.
%       This function uses the MIMO_MaximumLikelihood toolbox.
%										
%       References:
%
%                   Pintelon, R., J. Schoukens, and P. Guillaume (2007). Box-Jenkins identification revisited - Part III: multivariable 
%                   systems, Automatica, vol. 43, no. 5, pp. 868-875.
%
%                   Pintelon, R., J. Schoukens, G. Vandersteen, and K. Barbé (2010). Estimation of nonparametric noise and FRF models 
%                   for multivariable systems - Part II: extensions, applications, Mechanical Systems and Signal Processing, vol. 24, 
%                   no. 3, pp. 596-616.
%
%                   Pintelon, R., and J. Schoukens (2012). System Identification: A Frequency Domain Approach, second edition, 
%                   IEEE Press-Wiley, Piscataway (USA). 
%
%
%	Output parameters
%
%		ThetaPlant			=	estimated value plant, noise, and initial conditions parameters
%								ThetaPlant = struct('A', [], 'B', [], 'Ig', [], 'C',[],'D',[], 'Ih', [])
%									ThetaPlant.A     =   1 x (OrderA+1)
%                                                           ThetaPlant.A(r) = coefficient a(r-1) of Omega^(r-1) 
%									ThetaPlant.B     =   ny x nu x (OrderB+1)
%                                                           ThetaPlant.B(i,j,r) = coefficient b(i,j,r-1) of Omega^(r-1)
%									ThetaPlant.Ig    =   ny x (OrderIg+1) 
%
%
%	Input parameters
%
%		data				=	structure containing the non-parametric data required for the identification
%									data.Y		=	output DFT spectrum ny x 1, dimensions ny x number of frequencies
%									data.U		=	known input DFT spectrum nu x 1, dimensions nu x number of frequencies
%									data.freq	=	vector of frequency values (Hz), dimension: number of frequencies x 1
%									data.Ts		=	sampling time (s)
%									data.Gc		=	controller transfer function (feedback dynamics), must be given in case of a 
%                                                   system operating in feedback. Zero or empty if unknown or not present. 
%                                                   (optional, default zero) 
%                                                       size nu x ny x number of frequencies
%
%		Sel					=	structure with fields 'A', 'B', 'Ig'
%									Sel = struct('A', [], 'B', [], 'Ig', [], 'C', [], 'D', [], 'Ih', [])
%									Sel.A = 1 x (OrderA+1)
%										Sel.A(r) = 1 if coeff. a(r-1) is unknown
%										Sel.A(r) = 0 if coeff. a(r-1) = 0
%									Sel.B = ny x nu x (OrderB+1)
%										Sel.B(i,j,r) = 1 if coeff. b(i,j,r-1) is unknown
%										Sel.B(i,j,r) = 0 if coeff. b(i,j,r-1) = 0
%									Sel.Ig = ny x (OrderIg+1)
%										Sel.Ig(i,r) = 1 if coeff. ig(i,r-1) is unknown
%										Sel.Ig(i,r) = 0 if coeff. ig(i,r-1) = 0
%
%		ModelVar			=	contains the information about the model to be identified
%								structure with fields 'Transient', 'ThePlane', 'TheModel', 'RecipPlant', 'RecipNoise'
%									ModelVar = struct('Transient', [], 'PlantPlane', [], 'Struct', [], 'RecipPlant',[], 'RecipNoise')
%									ModelVar.Transient		=	1 then the initial conditions of the plant and/or noise are estimated
%									ModelVar.PlantPlane		=	plane of the plant model
%																	's':	continuous-time;
%																	'w':	sqrt(s)-domain
%																	'z':	discrete-time;
%									ModelVar.RecipPlant		=	1 if plant model is reciprocal: G(i,j) = G(j,i)
%
%		IterVar				=	contains the information about the minimization procedure
%								structure with fields 'LM', 'MaxIter', 'TolParam', 'TolCost', 'TraceOn'
%									IterVar = struct('LM', [], 'MaxIter', [], 'TolParam', [], 'TolCost', [], 'TraceOn', [])
%									IterVar.LM 			=	1 then the Levenberg-Marquardt minimization scheme is used
%									IterVar.MaxIter 	=	maximum number of itterations of the minimization procedure
%									IterVar.TolParam 	=	relative precision on parameter estimates
%									IterVar.TolCost 	=	relative precision on cost function
%									IterVar.TraceOn 	=	1 then output iterations (optional)
%
%		FigNum				=	number figure if a plot must be shown
%                               (optional, default 0: no plot is shown)
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, 20 October 2011 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add number of inputs and outputs to ModelVar
ModelVar.ny = size(data.Y,1);
ModelVar.nu = size(data.U,1);
nu = ModelVar.nu;
ny = ModelVar.ny;
% number of frequencies
F = length(data.freq);
% parametric estimate transient term
ModelVar.Transient = 1;
% data.freq should be a row vector
data.freq = data.freq(:).';

% feedback dynamics
try 
    if isempty(data.Gc)
        Feedback = 0;       % the plant operates in open loop
    else
        Feedback = 1;       % the plant operates in closed loop
    end % if
catch
    Feedback = 0;           % the plant operates in open loop
end % try

% calculation reference signal in case of closed loop operation
if Feedback == 1   
    dummy = zeros(1, ny, F);
    dummy(1,:,:) = data.Y;
    GcY = squeeze(sum(data.Gc .* repmat(dummy, [nu, 1, 1]), 2));
    if nu == 1
        GcY = GcY.';
    end % if single input
    data.R = GcY + data.U;
    ModelVar.Struct = 'EIV';   
else % plant operates in open loop  
    data.R = [];
    ModelVar.Struct = 'OE';    
end % if plant operates in closed loop

% figure number
try 
    if isempty(FigNum)
        FigNum = 0;
    end % if
catch
    FigNum = 0;
end % try


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation sample noise covariances of the input-output DFT spectra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% struct needed by the LocalPolyAnal routine:
% 2nd order local polynomial approximation; 6 degrees of freedom;
% transient elimination on; and all frequencies are used
method = struct('order', 2, 'dof', 6, 'transient', 1, 'step', 1); 

switch Feedback
    
    case 0  % open loop
        
        [CY, Y, TY, G, CvecG, dof, CL] = LocalPolyAnal(data, method);
        
        % data for parametric estimation: the original input-output
        % DFT spectra are used in order to have a correct parametric 
        % estimate of the plant transient term
        data.CY = CY.n;                     % sample noise covariance
        data.CU = zeros(nu, nu, F);
        data.CYU = zeros(ny, nu, F);
        
    case 1  % closed loop

        dataEIV = data;
        dataEIV.U = data.R;
        dataEIV.Y = [data.Y; data.U];

        % generalized sample mean and sample covariance Z = [Y; U]
        % and estimate FRF from R to Z and its covariance
        [CZ, Z, TZ, GRZ, CvecGRZ, dof, CL] = LocalPolyAnal(dataEIV, method); 
        
        % FRM from U to Y and its covariance
        [G, CvecG] = FRF_EIV(GRZ, CvecGRZ);

        % data for parametric estimation: the original input-output
        % DFT spectra are used in order to have a correct parametric 
        % estimate of the plant transient term
        data.CY = CZ.n(1:ny,1:ny,:);
        data.CU = CZ.n(ny+1:end,ny+1:end,:);
        data.CYU = CZ.n(1:ny,ny+1:end,:);
    
end % switch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation parametric plant model % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% starting values plant model: GTLS estimate followed by BTLS estimate
ThetaWGTLS = MIMO_WGTLS(data, Sel, ModelVar);
ThetaBTLS = MIMO_BTLS(data, Sel, ThetaWGTLS, ModelVar, IterVar);

% SML estimate
ThetaPlant = MIMO_ML(data, Sel, ThetaBTLS, ModelVar, IterVar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison parametric and non-parametric estimates plant FRM % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FigNum > 0
       
    % powers of jw, sqrt(jw), or exp(-jw*Ts)
    switch ModelVar.PlantPlane
        case 'z', x.Plant = exp(-sqrt(-1)*2*pi*data.freq.'*data.Ts);
        case 's', x.Plant = sqrt(-1)*2*pi*data.freq.';
        case 'w', x.Plant = (sqrt(-1)*2*pi*data.freq.').^0.5;
    end

    % estimated plant transfer function
    PolyTransPlant = MIMO_ML_CalcPolyTrans(ThetaPlant, x);

    % estimated variances FRM entries: keep the diagonal elements CvecG only
    F = length(data.freq);
    varG = zeros(ny, nu, F);
    for kk=1:F
        varG(:, :, kk) = reshape(diag(CvecG(:, :, kk)), [ny, nu]);
    end % kk
    
    % comparison parametric and nonparametric plant FRMs
    figure(FigNum);
    mm = 0;
    for jj = 1:ny
        for ii = 1:nu
            mm = mm+1;
            subplot(ny, nu, mm)
            plot(data.freq, db(squeeze(G(jj,ii,:))), 'k', data.freq, db(squeeze(PolyTransPlant.G(jj,ii,:))), 'r', ...
                 data.freq, db(squeeze(G(jj,ii,:)-PolyTransPlant.G(jj,ii,:))), 'k--', ...
                 data.freq, db(squeeze(varG(jj,ii,:)))/2, 'r--');
            xlabel('Frequency (Hz)')
            ylabel('Amplitude (dB)')
        end % ii
    end % jj
    subplot(ny, nu, ceil(nu/2))
    title('data: black; initial model: red; diff.: black --; var(data): red --');
    zoom on; shg

end % if plot figure