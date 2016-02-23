function [M] = mergeband(varargin)
%MERGEBAND - merge multiple frequency bands data.
%       M = MERGEBAND(B1,B2,B3,...,nrofb)
% Bx    : structures to concatenate
% nrofb : number of transient bands to remove
% M     : structure with merged time & freq data
% note  : time data are concatenated in cells
%         freq data are concatenated in vectors
% author: Thomas Beauduin, University of Tokyo, 2016

if isvector(varargin{1}), input=varargin{1}; end
if iscell(varargin{1}), input=cell2mat(varargin{1}); end
for k=1:nargin-1
    if isstruct(varargin{k}), input(k)=varargin{k}; end
    if iscell(input(k)), input(k)=cell2mat(input(k)); end
end

% PRETREAT
% pretreat data vectors of given structures
nrofb = varargin{end};
for k=1:length(input)
    B = input(k);
    fn = fieldnames(B);
    for m=1:numel(fn)
       if length(B.(fn{m}))==B.nroft
           if ~strcmp(fn{m},'time')
               [B.(fn{m}),B.time]=pretreat(B.(fn{m}),B.nrofs,B.fs,nrofb);
           end
       end
    end
end

% TIME2FRF
% convert data to frequency domain and merge bands
i=1;
y{i}='pos_t';   x{i}='iq_ad';  n{i}='pt'; i=i+1;
y{i}='theta_m'; x{i}='iq_ad';  n{i}='pm'; i=i+1;
y{i}='theta_s'; x{i}='iq_ad';  n{i}='ps'; i=i+1;
y{i}='disp_s2'; x{i}='iq_ad';  n{i}='ds'; i=i+1;
y{i}='acc_mx';  x{i}='iq_ad';  n{i}='am'; i=i+1;
y{i}='acc_sx';  x{i}='iq_ad';  n{i}='as'; i=i+1;
y{i}='acc_tx';  x{i}='iq_ad';  n{i}='at'; i=i+1;
y{i}='acc_tz';  x{i}='acc_tx'; n{i}='az'; i=i+1;

for k=1:length(input)
    B = input(k);
    for i=1:length(y)
        if isfield(B,y{i})
            [D.(strcat('X',n{i})){k,1},D.(strcat('Y',n{i})){k,1},...
             D.(strcat('FRFs',n{i})){k,1},D.(strcat('FRFn',n{i})){k,1},...
             D.('freq'){k,1},D.(strcat('sX2',n{i})){k,1},...
             D.(strcat('sY2',n{i})){k,1},D.(strcat('cXY',n{i})){k,1},...
             D.(strcat('sCR',n{i})){k,1}] ...
            = time2frf_ml(B.(x{i}),B.(y{i}),B.fs,B.fl,B.fh,B.nrofs); 
        end
    end
end

fn=fieldnames(D);
for m=1:numel(fn)
    if iscell(D.(fn{m})), M.(fn{m})=cell2mat(D.(fn{m}));
    else                  M.(fn{m})=D.(fn{m});
    end
end

end

