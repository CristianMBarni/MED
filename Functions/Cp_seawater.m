function [Cp] = Cp_seawater(T,X)

% T is temperature in ºC
% range between 20-180ºC
% X is salinity in ppm

% s is the salinity in mg/kg (same as ppm)

s = X;

A = 4206.8 -6.6197*s +1.2288e-2*s^2;
B = -1.1262 +5.4178e-2*s -2.2719e-4*s^2;
C = 1.2026e-2 -5.3566e-4*s +1.8906e-6*s^2;
D = 6.8777e-7 +1.517e-6*s -4.4268e-9*s^2; 

Cp = (A+ B*T + C*T^2 + D*T^3)*10^(-3);

end

