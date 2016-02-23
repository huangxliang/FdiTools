function [Y] = wrap(X)
%WRAP - wrapping of imaginary phase in [-2pi,pi] band.
%
% X,Y   :
% author:
%

for i=1:length(X)
    if X(i) < -180, X(i:end)=X(i:end)+360; end
    if X(i) >  180, X(i:end)=X(i:end)-360; end
end
Y=X;

% 1. SHIFT REMOVAL
% remove wrapped data up-down shifts for visibility
shift = 160; j=0; 
first = 0; last=0;
for i=2:length(Y)
    jump_i=0; jump_j=0; end_j=0;
    if abs(Y(i)-Y(i-1))>shift
        if j~=0
            if Y(i)>Y(i-1), jump_i=+1;  % jump up
            else            jump_i=-1;
            end
        else                bgn_j=+1;   % begin of array
        end

        for j=i:length(Y)-1
            if abs(Y(j)-Y(j+1))>shift, break; end
        end
        
        if j~=length(Y)-1
            if Y(j)>Y(j+1), jump_j=-1;  % jump down
            else            jump_j=+1;
            end
        else                end_j=+1;   % end of array
        end
        
        if bgn_j == 1 && first == 0
            if Y(i)<Y(i-1), Y(1:i-1)=Y(1:i-1)-abs(Y(i)-Y(i-1));
            else            jump_i=+1;
            end
            first = 1;
        end
        if jump_i == -jump_j
            Yij = mean(Y(i:j));
            Y(i:j)=Y(i:j)-abs(Yij-Y(i-1));%mean([abs(Yij-Y(i-1)),abs(Yij-Y(j+1))]);
        end
        if end_j == 1 && last == 0
           Y(i:j+1)=Y(i:j+1)-abs(Y(i)-Y(i-1)); 
           last = 1;
        end
    end
end
end