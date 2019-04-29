function h = swenthalpy(T,P,xs)

Po = 100000;
ws = xs/1000; %kg/kg

% Fitting polynomials through EES for Temperatures between 5-95ºC
hw = 0.117603379 + 4.20689825*T - 0.000800678287*T^2 + 0.00000895478984*T^3 - 2.61832765E-08*T^4;           %[kJ/kg]
rhow = 999.95866 + 0.0426118305*T - 0.00725068716*T^2 + 0.0000385494141*T^3 - 1.19726675E-07*T^4;           %[kg/m^3]

%% CoolProp
% CoolProp regular use is too slow

% Tk = T + 273.153;
% hw = CoolProp.PropsSI('H', 'P', Po, 'T', Tk,'Water')/1000;
% rhow = CoolProp.PropsSI('D', 'P', Po, 'T', Tk, 'Water');

% Maybe through low-level interface it gets better
% disp([num2str('*********** TABULAR BACKENDS *****************')]);
% TAB = CoolProp.AbstractState.factory('TTSE&HEOS', 'Water');
% TAB.update(CoolProp.PT_INPUTS, Po, Tk);
% rhow = TAB.rhomass();
% hw = TAB.hmass()/1000;

%%
b1 = -2.348*10^4;
b2 = 3.152*10^5;
b3 = 2.803*10^6;
b4 = -1.446*10^7;
b5 = 7.826*10^3;
b6 = -4.417*10^1;
b7 = 2.139*10^(-1);
b8 = -1.991*10^4;
b9 = 2.778*10^4;
b10 = 9.728*10^1;
a1 = 8.020*10^2;
a2 = -2.001;
a3 = 1.677*10^(-2);
a4 = -3.060*10^(-5);
a5 = -1.613*10^(-5);

rhosw = rhow + ws*(a1) + (a2)*T + (a3)*T^2 + (a4)*T^3 + (a5)*T^2*ws;
vsw = 1/rhosw;
hswo = hw - ws*(b1 + b2*ws + b3*ws^2 + b4*ws^3 + b5*T + b6*T^2 + b7*T^3 + b8*ws*T + b9*ws^2*T + b10*ws*T^2)/1000;

h = hswo + vsw*(P - Po)/1000000;

end
