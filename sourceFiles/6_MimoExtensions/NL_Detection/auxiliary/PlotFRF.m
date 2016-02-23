function FigNum = PlotFRF(G, Gc, Grat, freq, Spacing, RealNum, FigNum);
%
%       FigNum = function PlotFRF(G, Gc, Grat, freq, Spacing, RealNum, FigNum);
%
%
%   OUTPUT
%
%       FigNum  =   number last plotted figure
%
%   INPUT
%
%       G       =   structure {'all', 'mean', 'stdn', 'stdNL'} containing mean value and standard deviations FRF at excited odd harmonics;
%                   G.all       =   Frequency Response Function (FRF) for all realisations and periods; size M x P x F
%                   G.mean      =   mean value FRF over the P consecutive periods; size M x F
%                   G.stdn      =   struct{'E', 'NE'} containing the noise standard deviation of the mean FRF value over the P consecutive periods
%                                   G.stdn.E    =   noise std calculated from the excited frequencies; size M x F
%                                   G.stdn.NE   =   noise std calculated from the non-excited odd frequencies (= extrapolation); size M x F
%                                   Note: a difference between G.stdn.E and G.stdn.NE indicates a non-stationary behaviour
%                   G.stdNL     =   total standard deviation mean value FRF over P consecutive periods:
%                                   G.stdNL.^2  =   noise variance + stochatic NL distortions; size M x F
%
%       Gc      =   structure {'all', 'mean', 'stdn', 'stdNL'} containing mean value and standard deviations of the corrected FRF at excited odd harmonics,
%                   Note: a first order correction for the time variation is applied to the FRF
%                   Gc.all      =   Frequency Response Function (FRF) for all realisations and periods; size M x P x F
%                   Gc.mean     =   mean value FRF over the P consecutive periods
%                   Gc.stdn     =   struct{'E', 'NE'} containing the noise standard deviation of the mean FRF value over the P consecutive periods
%                                   Gc.stdn.E    =   noise std calculated from the excited frequencies
%                                   Gc.stdn.NE   =   noise std calculated from the non-excited odd frequencies (= extrapolation)
%                                   Note: a difference between G.stdn.E and G.stdn.NE indicates a non-stationary behaviour
%                   Gc.stdNL    =   total standard deviation mean value FRF over P consecutive periods:
%                                   Gc.stdNL.^2  =   noise variance + stochatic NL distortions; size M x F
%
%       Grat    =   struct{'all', 'mean'} containg the ratios of the FRF's over consecutive periods
%                   Grat.all    =   ratio Gall(:, ii+1, :) ./ Gall(:, ii, :)
%                   Grat.mean   =   mean value Grat over consecutive periods
%
%       freq    =   struct {'E', 'NE'} containing the excited and non-excited frequencies
%                   freq.E  =   vector containing the excited odd harmonic frequencies
%                   freq.NE =   struct{'even', 'odd'}
%                               freq.NE.even    =   struct{'all', 'inband', outband'}
%                                                   freq.NE.even.all        =   all non-excited even harmonics
%                                                   freq.NE.even.inband     =   non-excited inband even harmonics
%                                                   freq.NE.even.outband    =   non-excited outband even harmonics
%                               freq.NE.odd    =   struct{'all', 'inband', outband'}
%                                                   freq.NE.odd.all         =   all non-excited odd harmonics
%                                                   freq.NE.odd.inband      =   non-excited inband odd harmonics
%                                                   freq.NE.odd.outband     =   non-excited outband odd harmonics
%
%       Spacing =   'lin' or 'linear':      linear frequency axis in plot; default value
%                   'log' or 'logarithmic': logarithmic frequency axis
%
%       RealNum =   number of realisation that is shown in the single realisation figures; default = 1;
%
%       FigNum  =   start of figure numbers; default = 1;
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, April 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%

switch nargin
    case 3
        Spacing = 'lin';
    case 4
        RealNum = 1;
        FigNum = 1;
    case 5
        FigNum = 1;
end % switch
Spacing = lower(Spacing);

[M, F] = size(G.mean);
if M == 1
    RealNum = 1;
end % if

% plot the ratio of the FRF's over consecutive periods for one realisation
figure(FigNum); close; figure(FigNum);
switch Spacing
    case {'lin', 'linear'}
		subplot(211);
        plot(freq.E, db(squeeze(Grat.all(RealNum, :, :))), '+', freq.E, db(squeeze(Grat.mean(RealNum, :))), 'k');
        title('RATIO FRF: magnitude (dB)')
        subplot(212)
		plot(freq.E, angle(squeeze(Grat.all(RealNum, :, :)))*180/pi, '+', freq.E, angle(Grat.mean(RealNum, :))*180/pi,'k');
        title('RATIO FRF: phase (°)')
    case {'log', 'logarithmic'}
		subplot(211);
		semilogx(freq.E, db(squeeze(Grat.all(RealNum, :, :))), '+', freq.E, db(Grat.mean(RealNum, :)), 'k');
        title('magnitude (dB)')
        subplot(212)
		semilogx(freq.E, angle(squeeze(Grat.all(RealNum, :, :)))*180/pi, '+', freq.E, angle(Grat.mean(RealNum, :))*180/pi,'k');
        title('phase (°)')
end
zoom on;
shg

% plot the FRF and its uncertainties for one realisation
FigNum = FigNum + 1;
figure(FigNum); close; figure(FigNum);
switch Spacing
    case {'lin', 'linear'}
		subplot(211);
        plot(freq.E, db(G.mean(RealNum, :)), 'k', freq.E, db(G.stdn.E(RealNum,:)), 'g', freq.E, db(G.stdn.NE(RealNum,:)), 'g+', ...
             freq.E, db(G.stdNL(RealNum,:)), 'r')
        title('UNCORR.FRF: black: FRF; green: noise var. (-: excited; +: non-excited); red: total var.')
        subplot(212)
        plot(freq.E, db(Gc.mean(RealNum, :)), 'k', freq.E, db(Gc.stdn.E(RealNum,:)), 'g', freq.E, db(Gc.stdn.NE(RealNum,:)), 'g+', ...
             freq.E, db(Gc.stdNL(RealNum,:)), 'r')
        title('CORR. FRF: black: FRF; green: noise var. (-: excited; +: non-excited); red: total var.')
    case {'log', 'logarithmic'}
		subplot(211);
        semilogx(freq.E, db(G.mean(RealNum, :)), 'k', freq.E, db(G.stdn.E(RealNum,:)), 'g', freq.E, db(G.stdn.NE(RealNum,:)), 'g+', ...
                 freq.E, db(G.stdNL(RealNum,:)), 'r')
        title('UNCORR. FRF: black: FRF; green: noise var. (-: excited; +: non-excited); red: total var.')
        subplot(212)
        semilogx(freq.E, db(Gc.mean(RealNum, :)), 'k', freq.E, db(Gc.stdn.E(RealNum,:)), 'g', freq.E, db(Gc.stdn.NE(RealNum,:)), 'g+', ...
                 freq.E, db(Gc.stdNL(RealNum,:)), 'r')
        title('CORR. FRF: black: FRF; green: noise var. (-: excited; +: non-excited); red: total var.')
end
zoom on;
shg

if M > 1
	% rms averaging over the different realisations
    stdGm = struct('E', [], 'NE', []);
	Gmm= mean(G.mean, 1);
	Gsm = mean(G.stdNL.^2, 1).^0.5;
	stdGm.E = mean(G.stdn.E.^2, 1).^0.5;
	stdGm.NE = mean(G.stdn.NE.^2, 1).^0.5;
    
    stdGcm = struct('E', [], 'NE', []);
	Gcmm= mean(Gc.mean, 1);
	Gcsm = mean(Gc.stdNL.^2, 1).^0.5;
	stdGcm.E = mean(Gc.stdn.E.^2, 1).^0.5;
	stdGcm.NE = mean(Gc.stdn.NE.^2, 1).^0.5;

	FigNum = FigNum + 1;
    figure(FigNum); close; figure(FigNum);
	switch Spacing
        case {'lin', 'linear'}
			subplot(211);
            plot(freq.E, db(Gmm), 'k', freq.E, db(stdGm.E), 'g', freq.E, db(stdGm.NE), 'g+', freq.E, db(Gsm), 'r')
            title('UNCORR. FRF: black: FRF; green: noise var. (-: excited; +: non-excited); red: total var.')
            subplot(212)
            plot(freq.E, db(Gcmm), 'k', freq.E, db(stdGcm.E), 'g', freq.E, db(stdGcm.NE), 'g+', freq.E, db(Gcsm), 'r')
            title('CORR. FRF: black: FRF; green: noise var. (-: excited; +: non-excited); red: total var.')
        case {'log', 'logarithmic'}
			subplot(211);
            semilogx(freq.E, db(Gmm), 'k', freq.E, db(stdGm.E), 'g', freq.E, db(stdGm.NE), 'g+', freq.E, db(Gsm), 'r')
            title('UNCORR. FRF: black: FRF; green: noise var. (-: excited; +: non-excited); red: total var.')
            subplot(212)
            semilogx(freq.E, db(Gcmm), 'k', freq.E, db(stdGcm.E), 'g', freq.E, db(stdGcm.NE), 'g+', freq.E, db(Gcsm), 'r')
            title('CORR. FRF: black: FRF; green: noise var. (-: excited; +: non-excited); red: total var.')
	end
	zoom on;
	shg
    
end % if