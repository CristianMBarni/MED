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

%% System Inputs
n = 6; % Number of effects

% Flags
CT = 7; % Cooling technology: 7-PF; 8-PF/TVC; 9-FF; 10-FF/TVC
TVC_Strategy = 2; % 1-Lowest Compression Ratio; 2-Ts_sat is user defined; 3-Highest CR
d_Cross_flow_flag = 1; % 0/1 - Cross flow (Parallel)
eph_flag = 0; % 0-No preheaters; 1-PH/NCG; 2-alternating PH/NCG; 3-PH; 4-alternating PH
n_ph_NCG = 0; % Effect which NCG will start to act
pre_PlateHTX_flag = 1; % 0/1- Plate heat exchangers of the feedwater entering the down condenser

% Actual inputs
Mm = 1; % Mass flow of motive steam
Q_Loss = 0.3; % Fraction of thermal losses in each effect
Ts_sat = 100; % Temperature of saturated motive steam
Tv1 = 65; % Vapor temperature in the first effect
Tvn = 40; % Vapor temperature in the last effect
Tf_dc_out = 35; % Temperature of feedwater leaving the down condenser
pre_Tsw_out = 37; % Temperature of the seawater leaving the plate HXT
Tsw = 25; % Temperature of seawater
Tv_Loss = 0.2; % Temperature losses by pressure losses (Heat transfer inside each effect is at Tv-Tv_loss)
Xsw = 0.042; % Salinity of seawater in weight %
Xb1 = 0.045; % Salinity in weight % to be assumed when starting iteration - Calculate mass flow rate of feedwater for the 1st effect
Xbn_max = 0.07; % Maximum allowed salinity of the last effect

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

% n = 12; % Number of effects
% Ts = 70; % ºC
% Md = 139; % kg/s
% Xf = 35000; % ppm or mg/kg
% X(n) = 72000; % ppm or mg/kg
% T(n) = 41; % ºC
% deltaT_loss = 0.2; % ºC
% Tf = 35; % ºC
% Tcw = 25; % ºC

%% Problem solving

% Some loop

    % Reset salinity of first effect
    Xb(1) = Xb1;

    % Set temperature & pressure profiles
    T(1:n) = something;
    P(1:n) = something;