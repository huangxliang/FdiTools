function Mat = Mat_Mult(Mat1, Mat2);
%
%       fast multiplication of two 3-dimensional arrays
%
%   Mat = Mat_Mult(Mat1, Mat2);
%
%
%   Output
%
%       Mat     =   matrix product of Mat1 and Mat2, size rows x columns x F array  
%
%
%   Input
%
%       Mat1    =   rows x pp x F array with F >> rows, pp 
%
%       Mat2    =   pp x columns x F array with F >> pp, columns
%
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows, pp, F] = size(Mat1);
columns = size(Mat2, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast matrix multiplication %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:rows
    
    % select iith row of Mat1
    Mat1ii = squeeze(Mat1(ii,:,:));
    if pp == 1
        Mat1ii = Mat1ii.';                                      % squeeze on 1 x 1 x F gives F x 1 !!!  
    end % if
    
    for jj = 1:columns
        
        % select the jjth column of Mat2
        Mat2jj = squeeze(Mat2(:,jj,:));
        if pp == 1
            Mat2jj = Mat2jj.';                                  % squeeze on 1 x 1 x F gives F x 1 !!!
        end % if
        Mat(ii,jj,:) = sum(Mat1ii.*Mat2jj, 1);
        
    end %jj    
end % ii
