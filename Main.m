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
Md = 1; % kg/s
Xf = 42000; % ppm or mg/kg
X(n) = 70000; % ppm or mg/kg
T(n) = 40; % ºC
deltaT_loss = 2; % ºC
Tf = 35; % ºC
Tcw = 25; % ºC

%% Problem solving
% Available values:
hs_vap = latent_heat_water_evaporation(Ts); % enthalpy/mass released by the motive steam

Tv(n) = T(n) - deltaT_loss; % Vapor temperature accounting for losses
hv_vap(n) = latent_heat_water_evaporation(Tv(n)); % enthalpy/mass necessary to condensate the vapor on the nth effect

B(n) = (Xf/(X(n)-Xf))*Md; % Mass of brine released in the last efect
Mf = Md + B(n); % Mass of feedwater necessary

deltaT_total = Ts-T(n); % Overall temperature difference

%% Heat transfer coefficients
U(1) = 2.4; % Initial Guess
for i = 2:n
    U(i) = 0.95*U(i-1);
end

%% Temperature profile
deltaT(1) = deltaT_total/(U(1)*sum(1./U));
T(1) = Ts - deltaT(1);
aux = T(n);
for i = 2:n
    deltaT(i) = deltaT(1)*U(1)/U(i);
    T(i) = T(i-1) - deltaT(i);
end

if abs(aux - T(n)) > 0.01
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

%% Latent Heat
for i = 1:n
    Tv(i) = T(i) - deltaT_loss;
    hv_vap(i) = latent_heat_water_evaporation(Tv(i));
end

%% Distillate flow rate
aux = 0;
for i = 1:n
    aux = aux + hv_vap(1)/hv_vap(i);
end
D(1) = Md/aux;

for i = 2:n
    D(i) = D(1)*hv_vap(1)/hv_vap(i);
end

%% Brine flow rate

B(1) = Mf - D(1);
aux = B(n);
for i = 2:n
    B(i) = B(i-1) - D(i);
end
if abs(aux - B(n)) > 0.01
    warning('Specified final brine flow rate is different from calculated')
end

%% Salt concentration profile
X(1) = Xf*Mf/B(1);
aux = X(n);
for i = 2:n
    X(i) = X(i-1)*B(i-1)/B(i);
end
if abs(aux - X(n)) > 0.01
    warning('Specified final salt concentration is different from calculated')
end

%% Areas
A(1) = D(1)*hv_vap(1)/(U(1)*(Ts-T(1)));
for i = 2:n
    A(i) = D(i)*hv_vap(i)/(U(i)*(deltaT(i)-deltaT_loss));
end

%% Convergence criteria
% while max(abs(A(1:end-1)-A(2:end))) > 0.0001
%     Am = mean(A);
%     for i = 1:n
%         deltaT(i) = deltaT(i)*A(i)/Am;
%     end
% end

