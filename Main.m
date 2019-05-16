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
% Ts = 100; % �C
% Md = 1; % kg/s
% Xf = 42000; % ppm or mg/kg
% X(n) = 70000; % ppm or mg/kg
% T(n) = 40; % �C
% deltaT_loss = 2; % �C
% Tf = 35; % �C
% Tcw = 25; % �C

n = 12; % Number of effects
Ts = 70; % �C
Md = 139; % kg/s
Xf = 35000; % ppm or mg/kg
X(n) = 72000; % ppm or mg/kg
T(n) = 41; % �C
deltaT_loss = 0.2; % �C
Tf = 35; % �C
Tcw = 25; % �C

A = design_step(n,Ts,Md,Xf,X,T,deltaT_loss,Tf,Tcw);