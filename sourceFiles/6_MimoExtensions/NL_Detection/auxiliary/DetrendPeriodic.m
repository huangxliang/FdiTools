function [y, TheTrend] = DetrendPeriodic(u, N);
%
%   function [y, TheTrend] = DetrendPeriodic(u, N);
%
%       linear detrending of a periodic signal in pieces of two consecutive periods
%
%
%   OUTPUT
%
%       y           =   detrended signal; size 1 x P * N
%       TheTrend    =   struct{'value', 'time'}
%                       TheTrend.value =   1 x P vector containing the trend
%                       TheTrend.time  =   1 x P vector containing the corresponding time in samples
%
%   INPUT
%
%       u       =   periodic signal with trend; size 1 x P * N
%       N       =   number of samples in one period
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, March 2006 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%

u = u(:);
P = length(u) / N;          % number of periods

% test if number of periods is even
if abs(P/2 - floor(P/2)) > eps
    EvenNumberOfPeriods = 0;
else
    EvenNumberOfPeriods = 1;
end % if

switch EvenNumberOfPeriods
    case 1
        u = reshape(u, [N, P]);                 % each column is a period
    case 0
         u = reshape(u, [N, P]);                % each column is a period  
         u(:, end:end+1) = u(:, end-1:end);     % the second last period is duplicated
         P = P + 1;                             % one period has been added
end % switch

% calculation trend
% assumption: linear in two consecutive periods
ytrend = mean(u, 1);                                                
TheSlope = diff(ytrend) / N;                                        % slope linear trends in consecutive periods
TheSlope = TheSlope(1:2:end);                                       % blocks of two consecutive periods
TheTime = ([0:1:2*N-1] - (N-1)/2).';                                % two periods
LinTrend = TheTime * TheSlope + repmat(ytrend(1:2:end), [2*N, 1]);  % linear trends to be subtracted from u

% reshape u in two consecutive periods
u = reshape(u, [2*N, P/2]);

% detrend
y = u - LinTrend;
y = reshape(y, [P*N, 1]);
y = y.';

if ~EvenNumberOfPeriods
    % remove duplicated second last period
    y(end-2*N+1:end-N) = [];
    ytrend(end-1) = [];
    P = P - 1;
end % if

TheTrend.value = ytrend;
TheTrend.time = [(N-1)/2: N: (N-1)/2 + (P-1)*N];
