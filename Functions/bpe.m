function [Tbpe] = bpe(T,X)

% T is temperature in ºC
% X is salt weight %

A = 8.325e-2 +1.883e-4*T +4.02e-6*T^2;
B = -7.625e-4 +9.02e-5*T -5.2e-7*T^2;
C = 1.522e-4 -3e-6*T -3e-8*T^2;

Tbpe = A*X + B*X^2+C*X^3;

end

