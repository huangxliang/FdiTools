function FigNum = Plot_Robust_NL_Anal(G, Y, U, freq, Spacing, FigNum);
%
%   function FigNum = Plot_Robust_NL_Anal(G, Y, U, freq, Spacing, FigNum);
%
%       Plot the result of the nonlinear analysis of the FRF and input/output spectra
%           - figures start with FigNum
%           - input parameters = output parameters of the function Robust_NL_Anal
%
%       Examples
%           1. FigNum = Plot_Robust_NL_Anal(G, Y, U, freq)
%           2. FigNum = Plot_Robust_NL_Anal(G, Y, U, freq, Spacing)
%           3. FigNum = Plot_Robust_NL_Anal(G, Y, U, freq, Spacing, FigNum)
%       
%
%   OUTPUT
%
%       FigNum  =   number last plotted figure
%
%   INPUT
%
%	    G		:	struct('mean', 'stdNL', 'stdn', 'stds')
%                   G.mean  :   mean value of FRF estimate over the M x P FRF measurements
%                   G.stdNL :   total standard deviation mean FRF (contribution noise + NL distortion)
%                   G.stdn  :   noise standard deviation mean FRF
%                   G.stds  :   standard deviation stochastic NL distortions w.r.t. one multisine experiment
%                               = order of magnitude of NL bias errors = the difference between the best linear approximation G
%                                 and the true underlying linear system
%
%       Y       :   struct('mean', 'abs', 'stdn', 'stdNL', 'stds'); can be empty 
%                   Y.mean	:   if Rall is available then
%                                    Y.mean is the mean value of the projected output Yall./Rall
%                                    over all periods and realisations
%                               else
%                                   Y.mean is empty except if M = 1 then Y.mean is the mean value of the
%                                   output over all periods
%                   Y.abs	:   if Rall is available then
%                                   Y.abs = abs(Y.mean) 
%                               else
%                                   rms value |Y(k)| over different realisations compensated for the noise variance 
%                   Y.stdNL :   total variance Y.mean (as complex number) 
%                   Y.stdn  :   noise variance Y.mean (as complex number)
%                   Y.stds  :   output stochastic nonlinear distortions w.r.t. one multisine realisation 
%
%       U       :   struct('mean', 'abs', 'stdn', 'stdNL', 'stds'); can be empty 
%                   U.mean	:   if Rall is available then
%                                    U.mean is the mean value of the projected output Uall./Rall
%                                    over all periods and realisations
%                               else
%                                   U.mean is empty except if M = 1 then U.mean is the mean value of the
%                                   output over all periods
%                   U.abs	:   if Rall is available then
%                                   U.abs = abs(U.mean) 
%                               else
%                                   rms value |U(k)| over different realisations compensated for the noise variance 
%                   U.stdNL :   total variance U.mean (as complex number) 
%                   U.stdn  :   noise variance U.mean (as complex number)
%                   U.stds  :   input stochastic nonlinear distortions w.r.t. one multisine realisation 
%
%       freq    =   vector containing the frequencies
%
%       Spacing =   'lin' or 'linear':      linear frequency axis in plot; default value
%                   'log' or 'logarithmic': logarithmic frequency axis
%
%       FigNum  =   start of figure numbers; default = 1;
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, December 2007 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%

switch nargin
    case 4
        Spacing = 'lin';
        FigNum = 1;
    case 5
        FigNum = 1;
end % switch
Spacing = lower(Spacing);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the input/output DFT spectra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(U)
    U = struct('mean', [], 'abs', [], 'stdn', [], 'stdNL', [], 'stds', []);    
end % if no input

if isempty(Y)
    Y = struct('mean', [], 'abs', [], 'stdn', [], 'stdNL', [], 'stds', []);    
end % if no input

if (~isempty(U.abs)) & (~isempty(Y.abs))
    
    figure(FigNum);
    FigNum = FigNum+1;
    
    subplot(211)        % input
    switch Spacing
        
        case {'lin', 'linear'}
            if isempty(U.stdNL)     % total std is not available
                plot(freq, db(U.abs), 'k+', freq, db(U.stdn), 'g');
            else                    % if total std is available
                plot(freq, db(U.abs), 'k+', freq, db(U.stdn), 'g', freq, db(U.stdNL), 'r', freq, db(U.stds), 'b');
            end % if
            
        case {'log', 'logarithmic'}
            if isempty(U.stdNL)     % total std is not available
                semilogx(freq, db(U.abs), 'k+', freq, db(U.stdn), 'g');
            else                    % if total std is available
                semilogx(freq, db(U.abs), 'k+', freq, db(U.stdn), 'g', freq, db(U.stdNL), 'r', freq, db(U.stds), 'b');
            end % if 
            
    end % switch
    title('INPUT SPECTRUM: +: signal; green: noise std; red: total std; blue: stoch. NL std');
 
    subplot(212)        % output
    switch Spacing
        
        case {'lin', 'linear'}
            if isempty(Y.stdNL)     % total std is not available
                plot(freq, db(Y.abs), 'k+', freq, db(Y.stdn), 'g');
            else                    % if total std is available
                plot(freq, db(Y.abs), 'k+', freq, db(Y.stdn), 'g', freq, db(Y.stdNL), 'r', freq, db(Y.stds), 'b');
            end % if
            
        case {'log', 'logarithmic'}
            if isempty(Y.stdNL)     % total std is not available
                semilogx(freq, db(Y.abs), 'k', freq, db(Y.stdn), 'g');
            else                    % if total std is available
                semilogx(freq, db(Y.abs), 'k', freq, db(Y.stdn), 'g', freq, db(Y.stdNL), 'r', freq, db(Y.stds), 'b');
            end % if 
            
    end % switch    
    title('OUTPUT SPECTRUM: +: signal; green: noise std; red: total std; blue: stoch. NL std');
    
end % non-empty input/output


%%%%%%%
% FRF %
%%%%%%%

figure(FigNum);

switch Spacing

    case {'lin', 'linear'}
        plot(freq, db(G.mean), 'k', freq, db(G.stdn), 'g', freq, db(G.stdNL), 'r', freq, db(G.stds), 'b');

    case {'log', 'logarithmic'}
        semilogx(freq, db(G.mean), 'k', freq, db(G.stdn), 'g', freq, db(G.stdNL), 'r', freq, db(G.stds), 'b');

end % switch
title('ROBUST black: FRF; green: noise std; red: total std; blue: stoch. NL std');
