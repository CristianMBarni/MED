function [P] = P_sat_water_vapor(T)

% T is temperature in ºC
% Range between 5-200ºC with 0.05% error from steam table values

Pc = 22089; % kPa
Tc = 647.286; % K

f(1) = -7.419242;
f(2) = 0.29721;
f(3) = -0.1155286;
f(4) = 0.008685635;
f(5) = 0.001094098;
f(6) = -0.00439993;
f(7) = 0.002520658;
f(8) = -0.000521868;

aux = 0;
for i = 1:8
    aux = aux + f(i)*(0.01*(T +273.15 -338.15))^(i-1);
end

P = Pc*exp((Tc/(T+273.15)-1)*aux);
end