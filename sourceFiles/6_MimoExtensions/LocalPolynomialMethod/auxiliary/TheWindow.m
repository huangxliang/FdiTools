function [y, scale] = TheWindow(N, WindowType, NumberOfPeriods);
%
%           Calculates a window of the form
%
%               y(t) = 1 + alpha * cos(2*pi*t/N) + beta * cos(6*pi*t/N)
%
%           with t = 0,1, ..., N-1. According to the choice of alpha and
%           beta one gets:
%               - the rectangular window: alpha = beta = 0
%               - the Hanning window (within a scale factor): alpha = -1, beta = 0 
%               - the minimum variance window: alpha = beta = -1/2 
%               - the minimum leakage window: alpha = -9/8, beta = 1/8
%
%           Notes:
%               - the mean value of the time domain signals should be removed 
%                 before applying the time domain window 
%               - at least two signal periods should be measured for all
%                 windows, except the rectangular window
%               - to recover the correct Fourier coefficients of periodic signals, the DFT
%                 should be divided by lambda1 = sum(y)
%               - to preserve the variance of random signals the DFT should
%                 be divided by lambda2 = sqrt(sum(y.^2))
%               - the SNR of measuring the Fourier coefficients is
%                 decreased by the factor sqrt(1 + alpha^2/2 + beta^2/2) = sqrt(N)*lambda2/lambda1
%
%   function [y, scale] = TheWindow(N, WindowType, NumberOfPeriods);
%
%
%   Output
%
%       y               =   window, size 1 x N
%       scale           =   struct('period', 'rand', 'snr')
%                               scale.period    :   scale factor DFT to recover the correct Fourier coefficients
%                                                   of periodic signals 
%                               scale.rand      :   scale factor DFT to recover the correct variance of random signals
%                               scale.snr       :   factor giving the decrease in signal-to-noise ratio of the measured 
%                                                   Fourier coeffients
%
%   Input
%
%       N               =   number of time domain samples
%       WindowType      =   'rect', 'poly'  :   rectangular window
%                           'hann'          :	hanning window
%                           'mvar           : 	minimum variance window
%                           'mleak'         :   minimum leakage window
%
%       NumberOfPeriods =   optional parameter for the 'mvar' and 'mleak' windows: number of periods within the window (default 2) 
%
% Rik Pintelon, 6 februari 2009 
% version 26 March 2009
%

% intitialisation variables
try
    if isempty(NumberOfPeriods)
        NumberOfPeriods = 2;
    end
catch
    NumberOfPeriods = 2;
end % try

if NumberOfPeriods < 2
    NumberOfPeriods = 2;
end % if

scale = struct('period', [], 'rand', [], 'snr', []);

switch WindowType
    
    case {'rect', 'poly'}
        alpha = 0;
        beta = 0;
        
    case 'hann'
        alpha = -1;
        beta = 0;
        
    case 'mvar'
        alpha = -1/2;
        beta = -1/2;
        
    case 'mleak'
        alpha = -9/8;
        beta = 1/8;
    
end % switch

if NumberOfPeriods == 2
    kk = 3;
else
    kk = 2;
end % if

t = [0:1:N-1];
y = 1 + alpha * cos(2*pi*t/N) + beta * cos(2*pi*kk*t/N);

scale.period = sum(y);
scale.rand = sqrt(sum(y.^2));
scale.snr = sqrt(N)*scale.rand/scale.period;    % = sqrt(1 + alpha^2 + beta^2)
