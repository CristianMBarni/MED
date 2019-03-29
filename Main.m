%% Initialization
clc
close all
clearvars
path(pathdef)
addpath('Functions')

%% Useful Information
% Available functions:
% bpe
% Cp_seawater
% enthalpy_saturated_liquid_water
% enthalpy_saturated_water_vapor
% latent_heat_water_evaporation
% nea
% P_sat_water_vapor
% T_sat_water_vapor
% uc
% ue

%% Initial variables setup
n = 6; % Number of effects
Ts = 100; % ºC
Tn = 40; % ºC
Xn = 70000; % ppm or mg/kg
Xf = 42000; % ppm or mg/kg
Md = 1; % kg/s
dT_loss = 2; % ºC
Tf = 35; % ºC
Tcw = 25; % ºC

%% Problem solving

% Available values:
hs_vap = latent_heat_water_evaporation(Ts); % enthalpy/mass released by the motive steam

Tv(n) = Tn - dT_loss; % Vapor temperature accounting for losses
hv_vap(n) = latent_heat_water_evaporation(Tv(n)); % enthalpy/mass necessary to condensate the vapor on the nth effect

B(n) = (Xf/(Xn-Xf))*Md; % Mass of brine released in the last efect
Mf = Md + B(n); % Mass of feedwater necessary

deltaT_total = Ts-Tn; % Overall temperature difference

%%
U(1) = 2.4; % Initial Guess
for i = 2:n
    U(i) = 0.95*U(i-1);
end

deltaT(1) = deltaT_total/(U(1)*sum(1./U));
T(1) = Ts - deltaT(1);
for i = 2:n
    deltaT(i) = deltaT(1)*U(1)/U(i);
    T(i) = T(i-1) - deltaT(i);
end

if T(n) - Tn > 0.01
    warning('Specified final temperature is different from calculated')
end

% It should be noted that the temperature drop per effect increases as the
% effect temperature is reduced, i.e., dT1 > dT2 > dT3 > dT4 > dT5 > dT6.
% This is dictated by:
% - Constant heat transfer area,
% - Lower overall heat transfer coefficients at lower temperatures, and
% - Constant thermal loads in all effects.
% Therefore, the increase of the temperature drop at lower temperatures
% compensates the decrease in the overall heat transfer coefficient. 

