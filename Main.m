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
% n = 6; % Number of effects
% Ts = 100; % ºC
% Md = 1; % kg/s
% Xf = 42000; % ppm or mg/kg
% X(n) = 70000; % ppm or mg/kg
% T(n) = 40; % ºC
% deltaT_loss = 2; % ºC
% Tf = 35; % ºC
% Tcw = 25; % ºC

n = 12; % Number of effects
Ts = 70; % ºC
Md = 139; % kg/s
Xf = 35000; % ppm or mg/kg
X(n) = 72000; % ppm or mg/kg
T(n) = 41; % ºC
deltaT_loss = 0.2; % ºC
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
% U(1) = 2.4; % Initial Guess
U(1) = ue(Ts);
for i = 2:n
    U(i) = 0.95*U(i-1);
end

%% Initial temperature profile
deltaT(1) = deltaT_total/(U(1)*sum(1./U));
for i = 2:n
    deltaT(i) = deltaT(1)*U(1)/U(i);
end

% It should be noted that the temperature drop per effect increases as the
% effect temperature is reduced, i.e., dT1 > dT2 > dT3 > dT4 > dT5 > dT6.
% This is dictated by:
% - Constant heat transfer area,
% - Lower overall heat transfer coefficients at lower temperatures, and
% - Constant thermal loads in all effects.
% Therefore, the increase of the temperature drop at lower temperatures
% compensates the decrease in the overall heat transfer coefficient.
[D,B,X,T,Tv,hv_vap,A] = MED_equations(n,Md,Mf,Xf,U,deltaT,deltaT_loss,Ts);


%% Convergence criteria
iteration = 1;
while max(abs(A(1:end-1)-A(2:end))) > 0.0001
    disp(['Iteration ' num2str(iteration)])
    Am = mean(A);
    for i = 1:n
        deltaT(i) = deltaT(i)*A(i)/Am;
    end
    
%     U(1) = ue(T(1));
%     for i = 2:n
%         U(i) = 0.95*U(i-1);
%     end
    
    [D,B,X,T,Tv,hv_vap,A] = MED_equations(n,Md,Mf,Xf,U,deltaT,deltaT_loss,Ts);
    iteration = iteration + 1;
end

Ms = D(1)*hv_vap(1)/hs_vap;
PR = Md/Ms;

Qc = D(n)*hv_vap(n);
LMTDc = (Tf-Tcw)/log((T(n)-deltaT_loss-Tcw)/(T(n)-deltaT_loss-Tf));
Uc = uc(Tv(n));
Ac = Qc/(Uc*LMTDc);

sA = (sum(A)+Ac)/Md;

Cp = 4.2;
Mcw = D(n)*hv_vap(n)/(Cp*(Tf-Tcw)) - Mf;