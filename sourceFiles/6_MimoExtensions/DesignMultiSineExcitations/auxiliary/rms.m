function y = rms(x, dim);
%
%   function y = rms(x, dim);
%
%   OUTPUT
%
%       y       =   rms value of x along the dimension dim
%
%
%   INPUT
%
%       x       =   matrix containing the signal
%       dim     =   Optional parameter: dimension of the matrix along which the rms value is calculated 
%                   Default: dim = 1
%
%   Rik Pintelon, October 2009
%

try
    if isempty(dim)
        dim = 1;
    end % empty Gc
catch
    dim = 1;
end % try

y = mean(abs(x.^2), dim).^0.5;
