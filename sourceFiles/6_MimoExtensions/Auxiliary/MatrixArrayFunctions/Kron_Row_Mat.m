function Mat = Kron_Row_Mat(Row1, Mat1);
%
%       fast Kronecker product of the transpose of a 3-dimensional vector and a 3-dimensional matrix 
%
%   Mat = Kron_Row_Mat(Row1, Mat1)
%
%
%   Output
%
%       Mat     =   Kronecker product product of Vec1.' and Mat1: Mat = kron(Row1, Mat1) 
%                   size rows x pp*columns x F array  
%
%
%   Input
%
%       Row1    =   1 x pp x F array with F >> pp 
%
%       Mat1    =   rows x columns x F array with F >> rows, columns
%
%
% Copyright (c) Rik Pintelon, Vrije Universiteit Brussel - dept. ELEC, November 2009
% All rights reserved.
% Software can be used freely for non-commercial applications only.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation of the variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pp = size(Row1, 2);
[rows, columns, F] = size(Mat1);
Mat = zeros(rows, pp*columns, F);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast Kronecker product %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the matrix Mat for all frequencies
for ii = 1:pp
    Vec1ii = Row1(1,ii,:);
    Sel_Columns = [(ii-1)*columns+1:ii*columns];
    Mat(:, Sel_Columns, :) = Mat1 .* repmat(Vec1ii, [rows, columns, 1]);
end % ii
