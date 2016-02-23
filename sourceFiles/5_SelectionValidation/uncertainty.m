function [ output_args ] = uncertainty( input_args )
% code from sam, has to be tested

% Uncertainty
l_I = zeros(numMeas,size(frf,2));
figure
for i=1:numMeas
    l_I(i,:) = abs(( mag_nls'-abs(frf(i,:)) ) ./ mag_nls');
    semilogx(freq,l_I(i,:))
    hold on
end
hold off


end

