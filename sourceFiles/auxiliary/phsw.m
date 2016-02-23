function [Y] = phsw(X)
%PHS phs of tranfer data

Y = unwrap(angle(X))*180/pi;

end

