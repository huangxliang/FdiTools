function DvecLambda = CalcDvecLambda(ny);
%
% function DvecLambda = CalcDvecLambda(ny);
%
%   Output parameter
%       DvecLambda  =   derivative vec(Lambda, where Lambda has dimensions ny x ny, w.r.t. its (ny^2 + ny)/2 independent entries,
%                       dimension ny^2 x (ny^2 + ny)/2.
%
%   Input parameter
%       ny          =   number of outputs
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, April 2005 
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%

% sequence independent entries Lambda(rr, ss)
ii = 0;
TheIndex = cell((ny^2 + ny)/2,1);
for cc = 1:ny
    for rr = cc:ny
        ii = ii + 1;
        TheIndex{ii} = [rr, cc];
    end % rr
end % cc

% % sequence independent imaginary parts Lambda(rr, cc)
% for cc = 1:ny
%     for rr = cc+1:ny
%         ii = ii + 1;
%         TheIndex{ii} = [rr, cc];
%     end % rr
% end % cc

% position Lambda(rr, cc) in vec(Lambda)
TheRow = zeros(ny, ny);
for cc = 1:ny
    for rr = 1:ny
        TheRow(rr, cc) = rr + ny*(cc-1);
    end % rr
end % cc

% derivative vec(Lambda) w.r.t its independent entries
iimax = (ny^2 + ny)/2;
DvecLambda =  zeros(ny^2, iimax);
for ii = 1:iimax
    j1 = TheRow(TheIndex{ii}(1,1), TheIndex{ii}(1,2));
    j2 = TheRow(TheIndex{ii}(1,2), TheIndex{ii}(1,1));
    DvecLambda(j1, ii) = 1;
    DvecLambda(j2, ii) = 1;
end % ii

% % derivative vec(Lambda) w.r.t its independent imaginary parts
% for ii = iimax+1:ny^2
%     j1 = TheRow(TheIndex{ii}(1,1), TheIndex{ii}(1,2));
%     j2 = TheRow(TheIndex{ii}(1,2), TheIndex{ii}(1,1));
%     DvecLambda(j1, ii) = sqrt(-1);
%     DvecLambda(j2, ii) = -sqrt(-1);
% end % ii
