function NonExcitedHarm = HarmonicContent(MeasHarm, ExcitedHarm);
%
%
%   function NonExcitedHarm = HarmonicContent(MeasHarm, ExcitedHarm);
%
%
%   OUTPUT
%
%       NonExcitedHarm      =   'odd' multisine: structure classfying the measured non-excited harmonics in {'even', 'odd'}
%                                   NonExcitedHarm.even     =   structure classifying the measured non-excited even harmonics in {'all', 'inband', 'outband'}
%                                                               NonExcitedHarm.even.all     =   all measured non-excited even harmonics
%                                                               NonExcitedHarm.even.inband  =   measured inband non-excited even harmonics
%                                                               NonExcitedHarm.even.outband =   measured outband non-excited even harmonics
%                                   NonExcitedHarm.odd      =   structure classifying the measured non-excited odd harmonics in {'all', 'inband', 'outband'}
%                                                               NonExcitedHarm.odd.all      =   all measured non-excited odd harmonics
%                                                               NonExcitedHarm.odd.inband   =   measured inband non-excited odd harmonics
%                                                               NonExcitedHarm.odd.outband  =   measured outband non-excited odd harmonics
%                               'full' multisine: structure classfying the measured non-excited harmonics in {'all', 'inband', 'outband'}
%                                   NonExcitedHarm.all      =   all measured non-excited harmonics 
%                                   NonExcitedHarm.inband   =   in-band non-excited harmonics 
%                                   NonExcitedHarm.outband  =   out-band non-excited harmonics 
%
%
%   INPUT
%
%       MeasHarm        =   all measured harmonics
%       ExcitedHarm     =   the excited odd harmonics; all even harmonics are non-excited
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, November 2005 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
% Version 5 December 2007
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation structures %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MeasHarm = MeasHarm(:);
ExcitedHarm = ExcitedHarm(:);

% test whether only the odd harmonics are excited
OddMultisine = isempty(find((ExcitedHarm - floor(ExcitedHarm/2)*2) == 0, 1));

if OddMultisine
    
    NonExcitedHarm = struct('even', [], 'odd', []);
    NonExcitedHarm.even = struct('all', [], 'inband', [], 'outband', []);
    NonExcitedHarm.odd = struct('all', [], 'inband', [], 'outband', []);
    
else % full multisine
    
    NonExcitedHarm = struct('all', [], 'inband', [], 'outband', []);    
    
end % if


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selection of the excited and non-excited harmonics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if OddMultisine
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % selection of all even and odd measured harmonics %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zz = MeasHarm(1)/2;
    if (zz - floor(zz)) == 0                            % first measured harmonic is even
        NonExcitedHarm.even.all = MeasHarm(1:2:end);
        AllOddMeasHarm = MeasHarm(2:2:end);
    else                                                % first measured harmonic is odd
        NonExcitedHarm.even.all = MeasHarm(2:2:end);
        AllOddMeasHarm = MeasHarm(1:2:end);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % selection of all non-excited odd harmonics %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NonExcitedHarm.odd.all = AllOddMeasHarm;
    F = length(ExcitedHarm);
    RemoveIndex = zeros(F,1);
    for ii = 1:F
        RemoveIndex(ii) = find(AllOddMeasHarm == ExcitedHarm(ii));
    end
    NonExcitedHarm.odd.all(RemoveIndex) = [];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % selection inband and outband odd non-excited harmonics %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NonExcitedHarm.odd.inband = NonExcitedHarm.odd.all;
    NonExcitedHarm.odd.outband = NonExcitedHarm.odd.all;
    SelectIndex = [find(NonExcitedHarm.odd.all < ExcitedHarm(1)); find(NonExcitedHarm.odd.all > ExcitedHarm(end))];
    NonExcitedHarm.odd.inband(SelectIndex) = [];
    NonExcitedHarm.odd.outband = NonExcitedHarm.odd.outband(SelectIndex);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % selection inband and outband non-excited even harmonics %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NonExcitedHarm.even.inband = NonExcitedHarm.even.all;
    NonExcitedHarm.even.outband = NonExcitedHarm.even.all;
    SelectIndex = [find(NonExcitedHarm.even.all < ExcitedHarm(1)); find(NonExcitedHarm.even.all > ExcitedHarm(end))];
    NonExcitedHarm.even.inband(SelectIndex) = [];
    NonExcitedHarm.even.outband = NonExcitedHarm.even.outband(SelectIndex);
    

else % full multisine
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % selection of all non-excited harmonics %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NonExcitedHarm.all = MeasHarm;
    F = length(ExcitedHarm);
    RemoveIndex = zeros(F,1);
    for ii = 1:F
        RemoveIndex(ii) = find(MeasHarm == ExcitedHarm(ii));
    end
    NonExcitedHarm.all(RemoveIndex) = [];    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % selection inband and outband non-excited harmonics %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NonExcitedHarm.inband = NonExcitedHarm.all;
    NonExcitedHarm.outband = NonExcitedHarm.all;
    SelectIndex = [find(NonExcitedHarm.all < ExcitedHarm(1)); find(NonExcitedHarm.all > ExcitedHarm(end))];
    NonExcitedHarm.inband(SelectIndex) = [];
    NonExcitedHarm.outband = NonExcitedHarm.outband(SelectIndex);
    
    
end % if OddMultisine
