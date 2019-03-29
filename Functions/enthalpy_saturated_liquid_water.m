function [Cp] = Cp_seawater(T,X)

% T is temperature in ºC
% range between 20-180ºC
% X is salinity in ppm



Cp = (A+ B*T + C*T^2 + D*T^3)*10^(-3);

end

