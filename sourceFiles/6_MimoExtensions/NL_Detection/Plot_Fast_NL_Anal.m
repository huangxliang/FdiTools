function FigNum = Plot_Fast_NL_Anal(Y, Yc, U, G, freq, Spacing, RealNum, FigNum);
%
%   function FigNum = Plot_Fast_NL_Anal(Y, Yc, U, G, freq, Spacing, RealNum, FigNum);
%
%       Plot the result of the nonlinear analysis of the input/output spectra
%           - figures start with FigNum
%           - the single realisation number RealNum
%           - if more than one realisation is available then the robust analysis
%             of the FRF is calculated and compared to the fast method mean
%             rms averaged over the different realisations
%       
%
%   OUTPUT
%
%       FigNum  =   number last plotted figure
%
%   INPUT
%
%       Y       =   struct{'mean', 'std'} containing mean value and standard deviation output spectrum
%                   Y.mean  =   struct{'E', 'NE'} containing the mean value of the output spectrum over P consecutive periods; size M x P x NumberOfFrequencies
%                               Y.mean.E    =   mean output spectrum at the excited (odd) harmonic frequencies
%                               Y.mean.NE   =   struct{'even', 'odd'} for ODD multisines
%                                               struct{'all', 'inband', outband'} for FULL multisines 
%                                               Y.mean.NE.even      =   struct {'all', 'inband', 'outband'}
%                                                                       Y.mean.NE.even.all      =   mean output spectrum at all the non-excited even harmonics
%                                                                       Y.mean.NE.even.inband   =   mean output spectrum at the inband non-excited even harmonics
%                                                                       Y.mean.NE.even.outband  =   mean output spectrum at the outband non-excited even harmonics
%                                               Y.mean.NE.odd       =   struct {'all', 'inband', 'outband'}
%                                                                       Y.mean.NE.odd.all       =   mean output spectrum at all the non-excited odd harmonics
%                                                                       Y.mean.NE.odd.inband    =   mean output spectrum at the inband non-excited odd harmonics
%                                                                       Y.mean.NE.odd.outband   =   mean output spectrum at the outband non-excited odd harmonics
%                                               Y.mean.NE.all       =   mean output spectrum at all non-excited harmonics 
%                                               Y.mean.NE.inband    =   mean output spectrum at the inband non-excited harmonics 
%                                               Y.mean.NE.outband	=   mean output spectrum at the outband non-excited harmonics 
%                   Y.std   =   struct{'E', 'NE'} containing the standard deviation of the mean value of the output spectrum over P consecutive periods
%                               Y.std.E     =   std output spectrum at the excited (odd) harmonic frequencies
%                               Y.std.NE    =   struct{'even', 'odd'} for ODD multisines 
%                                               struct{'all', 'inband', outband'} for FULL multisines 
%                                               Y.std.NE.even       =   struct {'all', 'inband', 'outband'}
%                                                                       Y.std.NE.even.all       =   std output spectrum at all the non-excited even harmonics
%                                                                       Y.std.NE.even.inband    =   std output spectrum at the inband non-excited even harmonics
%                                                                       Y.std.NE.even.outband   =   std output spectrum at the outband non-excited even harmonics
%                                               Y.std.NE.odd        =   struct {'all', 'inband', 'outband'}
%                                                                       Y.std.NE.odd.all        =   std output spectrum at all the non-excited odd harmonics
%                                                                       Y.std.NE.odd.inband     =   std output spectrum at the inband non-excited odd harmonics
%                                                                       Y.std.NE.odd.outband    =   std output spectrum at the outband non-excited odd harmonics
%                                               Y.std.NE.all        =   std output spectrum at all non-excited harmonics 
%                                               Y.std.NE.inband     =   std output spectrum at the inband non-excited harmonics 
%                                               Y.std.NE.outband	=   std output spectrum at the outband non-excited harmonics 
%
%       Yc      =   structure {'mean', 'std'} containing the inband mean value and standard deviation corrected output spectrum; size M x P x NumberOfFrequencies
%                   corrected for the presence of energy at non-excited input harmonics
%                   Yc.mean =   struct{'E', 'NE'} containing the inband mean corrected output spectrum
%                               Yc.mean.E   =   mean corrected output spectrum at the excited (odd) harmonics (same as uncorrected)
%                               Yc.mean.NE  =   struct{'even', 'odd'} for ODD multisines 
%                                               Yc.mean.NE.even =   mean corrected output spectrum at non-excited inband even harmonics
%                                               Yc.mean.NE.odd  =   mean corrected output spectrum at non-excited inband odd harmonics
%                                               mean corrected output spectrum at the non-excited harmonics for FULL multisines 
%                   Yc.std  =   struct{'E', 'NE'} containing the inband standard deviation of the corrected output spectrum
%                               Yc.std.E    =   std corrected output spectrum at excited odd harmonics (same as uncorrected)
%                               Yc.std.NE   =   struct{'even', 'odd'} for ODD multisines 
%                                               Yc.std.NE.even  =   std corrected output spectrum at non-excited inband even harmonics
%                                               Yc.std.NE.odd   =   std corrected output spectrum at non-excited inband odd harmonics
%                                               std corrected output spectrum at the non-excited harmonics for FULL multisines 
%
%       U       =   structure {'mean', 'std'} containing mean value and standard deviation input spectrum; size M x P x F
%                   U.mean  =   struct{'E', 'NE'} containing the mean value of the input spectrum over P consecutive periods
%                               U.mean.E    =   mean input spectrum at excited odd harmonic frequencies
%                               U.mean.NE   =   struct{'even', 'odd'} for ODD multisines 
%                                               struct{'all', 'inband', outband'} for FULL multisines 
%                                               U.mean.NE.even      =   struct {'all', 'inband', 'outband'}
%                                                                       U.mean.NE.even.all      =   mean input spectrum at all the non-excited even harmonics
%                                                                       U.mean.NE.even.inband   =   mean input spectrum at the inband non-excited even harmonics
%                                                                       U.mean.NE.even.outband  =   mean input spectrum at the outband non-excited even harmonics
%                                               U.mean.NE.odd       =   struct {'all', 'inband', 'outband'}
%                                                                       U.mean.NE.odd.all       =   mean input spectrum at all the non-excited odd harmonics
%                                                                       U.mean.NE.odd.inband    =   mean input spectrum at the inband non-excited odd harmonics
%                                                                       U.mean.NE.odd.outband   =   mean input spectrum at the outband non-excited odd harmonics
%                                               U.mean.NE.all       =   mean input spectrum at all non-excited harmonics 
%                                               U.mean.NE.inband    =   mean input spectrum at the inband non-excited harmonics 
%                                               U.mean.NE.outband	=   mean input spectrum at the outband non-excited harmonics 
%                   U.std   =   struct{'E', 'NE'} containing the standard deviation of the mean value of the input spectrum over P consecutive periods
%                               U.std.E     =   std input spectrum at excited odd harmonic frequencies
%                               U.std.NE    =   struct{'even', 'odd'} for ODD multisines 
%                                               struct{'all', 'inband', outband'} for FULL multisines 
%                                               U.std.NE.even       =   struct {'all', 'inband', 'outband'}
%                                                                       U.std.NE.even.all       =   std input spectrum at all the non-excited even harmonics
%                                                                       U.std.NE.even.inband    =   std input spectrum at the inband non-excited even harmonics
%                                                                       U.std.NE.even.outband   =   std input spectrum at the outband non-excited even harmonics
%                                               U.std.NE.odd        =   struct {'all', 'inband', 'outband'}
%                                                                       U.std.NE.odd.all        =   std input spectrum at all the non-excited odd harmonics
%                                                                       U.std.NE.odd.inband     =   std input spectrum at the inband non-excited odd harmonics
%                                                                       U.std.NE.odd.outband    =   std input spectrum at the outband non-excited odd harmonics
%                                               U.std.NE.all        =   std input spectrum at all non-excited harmonics 
%                                               U.std.NE.inband     =   std input spectrum at the inband non-excited harmonics 
%                                               U.std.NE.outband	=   std input spectrum at the outband non-excited harmonics 
%
%       G       =   structure {'all', 'mean', 'stdn', 'stdNL', 'stds'} containing mean value and standard deviations FRF at excited odd harmonics
%                   G.all   =   Frequency Response Function (FRF) for all realisations and periods; size M x P x F
%                   G.mean  =   mean value FRF over the P consecutive periods
%                   G.stdn  =   struct{'E', 'NE'} containing the noise standard deviation of the mean FRF value over the P consecutive periods
%                               G.stdn.E    =   noise std calculated from the excited frequencies
%                               G.stdn.NE   =   noise std calculated from the non-excited odd frequencies (= extrapolation)
%                               Note: a difference between G.stdn.E and G.stdn.NE indicates a non-stationary behaviour
%                   G.stdNL =   total standard deviation mean value FRF over P consecutive periods:
%                               G.stdNL.^2  =   noise variance + stochatic NL distortions
%                   G.stds =    standard deviation stochastic NL distortions
%                               G.stds  =   sqrt(|G.stdNL.^2 - G.stdn.^2|)
%                               Note: G.stds is a measure of the order of magnitude of the bias contribution of the NL distortions to the
%                                     best linear approximation
%
%       freq    =   struct {'E', 'NE'} containing the excited and non-excited frequencies
%                   freq.E  =   vector containing the excited (odd) harmonic frequencies
%                   freq.NE =   struct{'even', 'odd'} for ODD multisines
%                               struct{'all', 'inband', outband'} for FULL multisines 
%                               freq.NE.even        =   struct{'all', 'inband', outband'}
%                                                       freq.NE.even.all        =   all non-excited even harmonics
%                                                       freq.NE.even.inband     =   inband non-excited even harmonics
%                                                       freq.NE.even.outband    =   outband non-excited even harmonics
%                               freq.NE.odd         =   struct{'all', 'inband', outband'}
%                                                       freq.NE.odd.all         =   all non-excited odd harmonics
%                                                       freq.NE.odd.inband      =   inband non-excited odd harmonics
%                                                       freq.NE.odd.outband     =   outband non-excited odd harmonics
%                               freq.NE.all         =   all non-excited harmonics 
%                               freq.NE.inband      =   inband non-excited harmonics 
%                               freq.NE.outband     =   outband non-excited harmonics 
%
%       Spacing =   'lin' or 'linear':      linear frequency axis in plot; default value
%                   'log' or 'logarithmic': logarithmic frequency axis
%
%       RealNum =   number of realisation that is shown in the single realisation figures; default = 1;
%
%       FigNum  =   start of figure numbers; default = 1;
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, March 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 5 December 2007
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 5
        Spacing = 'lin';
    case 6
        RealNum = 1;
        FigNum = 1;
    case 7
        FigNum = 1;
end % switch
Spacing = lower(Spacing);

[M, F] = size(Y.mean.E);
if M == 1
    RealNum = 1;
end % if

% check if multisine is odd or full
% if Y.mean.NE.all exists then full multisine; otherwise odd multisine
try
    OddMultisine = isempty(Y.mean.NE.all);
catch
    OddMultisine = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the results of the robust method             %
% - noise variances of the mean value               %
% - stochastic NL distortion w.r.t. one realisation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if M > 1
    % robust estimate of the level of the NL distortion w.r.t one realisation = M * (varGmmNL - varGmmn)
	Gmm= mean(G.mean, 1);
	stdGmmNL = std(G.mean, 0, 1)/sqrt(M);
	stdGmmn = (mean(G.stdn.E.^2, 1) / M) .^ 0.5;
    stdGmms = (abs(stdGmmNL.^2 - stdGmmn.^2) * M).^0.5;
	
    figure(FigNum); close; figure(FigNum);
    switch Spacing
        case {'lin', 'linear'}
			plot(freq.E, db(Gmm),'k', freq.E, db(stdGmmNL),'r', freq.E, db(stdGmmn),'g', freq.E, db(stdGmms), 'b')
        case {'log', 'logarithmic'}
			semilogx(freq.E, db(Gmm),'k', freq.E, db(stdGmmNL),'r', freq.E, db(stdGmmn),'g', freq.E, db(stdGmms), 'b')
    end
	title('ROBUST black: FRF; green: noise var. mean; red: total var. mean; blue: NL dist.')
	zoom on;
	shg
else
    FigNum = FigNum-1;
end % if


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the measured input-output spectra of one multisine measurement %
% Notes:                                                              %
%   - number of the plotted realisation is defined by RealNum         %
%   - even NL distortions fall outside the band of interest           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigNum = FigNum + 1;
figure(FigNum); close; figure(FigNum);

subplot(211)
switch Spacing
    case {'lin', 'linear'}
        if OddMultisine
            plot(freq.E, db(U.mean.E(RealNum, :)), 'k+', freq.E, db(U.std.E(RealNum, :)), 'k', ...
                 freq.NE.odd.all, db(U.mean.NE.odd.all(RealNum, :)), 'ro', freq.NE.odd.all, db(U.std.NE.odd.all(RealNum, :)), 'r', ...
                 freq.NE.even.all, db(U.mean.NE.even.all(RealNum, :)), 'b*', freq.NE.even.all, db(U.std.NE.even.all(RealNum, :)), 'b');
        else % full multisine
            plot(freq.E, db(U.mean.E(RealNum, :)), 'k+', freq.E, db(U.std.E(RealNum, :)), 'k', ...
                 freq.NE.all, db(U.mean.NE.all(RealNum, :)), 'r^', freq.NE.all, db(U.std.NE.all(RealNum, :)), 'r');
        end % if OddMultisine
    case {'log', 'logarithmic'}
        if OddMultisine
            semilogx(freq.E, db(U.mean.E(RealNum, :)), 'k+', freq.E, db(U.std.E(RealNum, :)), 'k', ...
                 freq.NE.odd.all, db(U.mean.NE.odd.all(RealNum, :)), 'ro', freq.NE.odd.all, db(U.std.NE.odd.all(RealNum, :)), 'r', ...
                 freq.NE.even.all, db(U.mean.NE.even.all(RealNum, :)), 'b*', freq.NE.even.all, db(U.std.NE.even.all(RealNum, :)), 'b');
        else % full multisine
            semilogx(freq.E, db(U.mean.E(RealNum, :)), 'k+', freq.E, db(U.std.E(RealNum, :)), 'k', ...
                 freq.NE.all, db(U.mean.NE.all(RealNum, :)), 'r^', freq.NE.all, db(U.std.NE.all(RealNum, :)), 'r');
        end % if OddMultisine
end
if OddMultisine
    title('INPUT SPECTRUM +: signal; 0: odd NL; *: even NL; -: std')
else % full multisine
    title('INPUT SPECTRUM +: signal; triangle: NL; -: std')
end % if OddMultisine

subplot(212)
switch Spacing
    case {'lin', 'linear'}
        if OddMultisine
            plot(freq.E, db(Y.mean.E(RealNum, :)), 'k+', freq.E, db(Y.std.E(RealNum, :)), 'k', ...
                 freq.NE.odd.all, db(Y.mean.NE.odd.all(RealNum, :)), 'ro', freq.NE.odd.all, db(Y.std.NE.odd.all(RealNum, :)), 'r', ...
                 freq.NE.even.all, db(Y.mean.NE.even.all(RealNum, :)), 'b*', freq.NE.even.all, db(Y.std.NE.even.all(RealNum, :)), 'b');
        else % full multisine
            plot(freq.E, db(Y.mean.E(RealNum, :)), 'k+', freq.E, db(Y.std.E(RealNum, :)), 'k', ...
                 freq.NE.all, db(Y.mean.NE.all(RealNum, :)), 'r^', freq.NE.all, db(Y.std.NE.all(RealNum, :)), 'r');
        end % if OddMultisine
    case {'log', 'logarithmic'}
        if OddMultisine
            semilogx(freq.E, db(Y.mean.E(RealNum, :)), 'k+', freq.E, db(Y.std.E(RealNum, :)), 'k', ...
                 freq.NE.odd.all, db(Y.mean.NE.odd.all(RealNum, :)), 'ro', freq.NE.odd.all, db(Y.std.NE.odd.all(RealNum, :)), 'r', ...
                 freq.NE.even.all, db(Y.mean.NE.even.all(RealNum, :)), 'b*', freq.NE.even.all, db(Y.std.NE.even.all(RealNum, :)), 'b');
        else % full multisine
            semilogx(freq.E, db(Y.mean.E(RealNum, :)), 'k+', freq.E, db(Y.std.E(RealNum, :)), 'k', ...
                 freq.NE.all, db(Y.mean.NE.all(RealNum, :)), 'r^', freq.NE.all, db(Y.std.NE.all(RealNum, :)), 'r');
        end % if OddMultisine
end
if OddMultisine
    title('OUTPUT SPECTRUM +: signal; 0: odd NL; *: even NL; -: std')
else % full multisine
    title('OUTPUT SPECTRUM +: signal; triangle: NL; -: std')
end % if OddMultisine

zoom on;
shg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the output spectrum corrected for the presence of NL distortions at the input %
% This correction is made at the non-excited inband harmonics                        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigNum = FigNum + 1;
figure(FigNum); close; figure(FigNum);

switch Spacing
    case {'lin', 'linear'}
        if OddMultisine
            plot(freq.E, db(Yc.mean.E(RealNum, :)), 'k+', freq.E, db(Yc.std.E(RealNum, :)), 'k', ...
                 freq.NE.odd.inband, db(Yc.mean.NE.odd(RealNum, :)), 'ro', freq.NE.odd.inband, db(Yc.std.NE.odd(RealNum, :)), 'r', ...
                 freq.NE.even.inband, db(Yc.mean.NE.even(RealNum, :)), 'b*', freq.NE.even.inband, db(Yc.std.NE.even(RealNum, :)), 'b');
        else % full multisine
            plot(freq.E, db(Yc.mean.E(RealNum, :)), 'k+', freq.E, db(Yc.std.E(RealNum, :)), 'k', ...
                 freq.NE.inband, db(Yc.mean.NE(RealNum, :)), 'r^', freq.NE.inband, db(Yc.std.NE(RealNum, :)), 'r');
        end % if OddMultisine
    case {'log', 'logarithmic'}
        if OddMultisine
            semilogx(freq.E, db(Yc.mean.E(RealNum, :)), 'k+', freq.E, db(Yc.std.E(RealNum, :)), 'k', ...
                 freq.NE.odd.inband, db(Yc.mean.NE.odd(RealNum, :)), 'ro', freq.NE.odd.inband, db(Yc.std.NE.odd(RealNum, :)), 'r', ...
                 freq.NE.even.inband, db(Yc.mean.NE.even(RealNum, :)), 'b*', freq.NE.even.inband, db(Yc.std.NE.even(RealNum, :)), 'b');
        else % if full multisine
            semilogx(freq.E, db(Yc.mean.E(RealNum, :)), 'k+', freq.E, db(Yc.std.E(RealNum, :)), 'k', ...
                 freq.NE.inband, db(Yc.mean.NE(RealNum, :)), 'r^', freq.NE.inband, db(Yc.std.NE(RealNum, :)), 'r');
        end % if OddMultisine
end
if OddMultisine
    title('CORRECTED OUTPUT SPECTRUM +: signal; 0: odd NL; *: even NL; -: std')
else % full multisine
    title('CORRECTED OUTPUT SPECTRUM +: signal; triangle: NL; -: std')
end % if OddMultisine

zoom on;
shg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the predicted distortion on FRF from the level of the non-excited (odd) harmonics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigNum = FigNum + 1;
figure(FigNum); close; figure(FigNum);
switch Spacing
    case {'lin', 'linear'}
        plot(freq.E, db(G.mean(RealNum, :)), 'k', freq.E, db(G.stdn.E(RealNum,:)), 'g', ...
             freq.E, db(G.stdNL(RealNum,:)), 'r')
    case {'log', 'logarithmic'}
        semilogx(freq.E, db(G.mean(RealNum, :)), 'k', freq.E, db(G.stdn.E(RealNum,:)), 'g', ...
                 freq.E, db(G.stdNL(RealNum,:)), 'r')
end
title('FAST - ONE REAL.  black: FRF; green: noise std; red: total std')
zoom on;
shg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the estimated noise-levels on the FRF measurement for one multisine %
% experiment obtained via                                                  % 
%   - excited harmonics                                                    % 
%   - non-excited harmonics                                                % 
% If both differ then either ther is a slip between the generator and      %
% acquisition unit, or the system is time-varying                          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigNum = FigNum + 1;
figure(FigNum); close; figure(FigNum);

switch Spacing
    case {'lin', 'linear'}
        plot(freq.E, db(G.mean(RealNum, :)), 'k', freq.E, db(G.stdn.E(RealNum,:)), 'g', ...
             freq.E, db(G.stdn.NE(RealNum,:)), 'g+')
    case {'log', 'logarithmic'}
        semilogx(freq.E, db(G.mean(RealNum, :)), 'k', freq.E, db(G.stdn.E(RealNum,:)), 'g', ...
                 freq.E, db(G.stdn.NE(RealNum,:)), 'g+')
end
title('FAST - ONE REAL.  black: FRF; green: noise var. using excited (-) and non-excited (+)')

zoom on;
shg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the fast and the robust method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if M > 1
	% rms averaging over the different realisations to compare information extracted by looking in between the lines
	% to the robust method; the total variance and noise variance of the fast method are rescaled to the mean over M realisations
	GNLm = mean(G.stdNL.^2 / M, 1).^0.5;
	stdGm = mean(G.stdn.E.^2 / M, 1).^0.5;
    Gsm = mean(G.stds.^2, 1).^0.5;
	
	FigNum = FigNum + 1;
    figure(FigNum); close; figure(FigNum);
	switch Spacing
        case {'lin', 'linear'}
        	plot(freq.E, db(Gmm), 'k', freq.E, db(stdGmmn), 'g', freq.E, db(stdGmmNL), 'r', ...
                 freq.E, db(stdGmms), 'b', freq.E, db(stdGm), 'g+', freq.E, db(GNLm), 'r+', freq.E, db(Gsm), 'b+')
        case {'log', 'logarithmic'}
        	semilogx(freq.E, db(Gmm), 'k', freq.E, db(stdGmmn), 'g', freq.E, db(stdGmmNL), 'r', ...
                     freq.E, db(stdGmms), 'b', freq.E, db(stdGm), 'g+', freq.E, db(GNLm), 'r+', freq.E, db(Gsm), 'b+')
	end
	title('COMP. ROBUST (-) and FAST (+) black: FRF; green: noise var; red: total var; blue: stoch. NL')
    zoom on;
    shg
end % if
