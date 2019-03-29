function [Tnea] = nea(Tb1,Tb2,X2)

% Tb1 is the boiling brine temperature on the previous effect
% Tb2 is the boiling brine temperature on the current effect
% Miyatake et al. (1973)

Tv = Tb2 - bpe(Tb2,X2);

Tnea = 33*(Tb1-Tb2)^0.55/Tv;

end

