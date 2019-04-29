function [H] = enthalpy_saturated_liquid_water(T)

% T is temperature in ºC
% range between 5-200ºC

H = -0.033635409 + 4.207557011.*T - 6.200339e-4.*T.^2 + 4.459374e-6.*T.^3;
end

