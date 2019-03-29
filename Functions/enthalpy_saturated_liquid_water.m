function [Cp] = Cp_seawater(T,X)

% T is temperature in �C
% range between 20-180�C
% X is salinity in ppm



Cp = (A+ B*T + C*T^2 + D*T^3)*10^(-3);

end

