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

%% Simplified MED initial variables setup - JUST FOR REFERENCE
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
Ms = 1; % Mass flow of motive steam
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


%% Problem solving
while(true)
    
    % Set/Reset 1st effect salinity
    Xb(1) = Xb1;
    
    %% Vapor and brine temperature profile across effects
    % Set temperature and pressure profile
    Tv(1) = Tv1;
    Tv(n) = Tvn;
    delta_Tv = (Tv(1) - Tv(n))/(n-1);
    Tv = Tv(1):-delta_Tv:Tv(n);
    
    Pv = P_sat_water_vapor(Tv);
    
    % Set temperature profile of the brine
    Tbpe = bpe(Tv,Xb(1));
    Tb = Tv + Tbpe;
    
    % ???? - Difference between Tb(i) and Tf(i) entering the effects - ????
    Delta_Tf_iph = Tb(n) - Tf_dc_out;
    
    % Assume HTX temp. inside the effect to be the average
    Delta_T_htx(2:n) = (Tb(2:n)+Tv(2:n))/2;
    
    % Calculate subcooled temperature of distillate collected from 1st effect
    Delta_T_htx(1) = sum(Delta_T_htx(2:n))/(n-1);
    Ts_sub = Tb(1) + Delta_T_htx(1);
    
    %%
    % Difference between Tb and Tf entering the effects
    Delta_Tf_iph = Tb(n) - Tf_dc_out;
    
    % Power available for 1st effect
    delta_Hs = latent_heat_water_evaporation(Ts_sat); % For saturated motive steam
    Q(1) = Ms*(delta_Hs);
    
    %% Temperature profile of feedwater input across effects
    Tf_eph(n) = Tf_dc_out;
    
    % There are no preheaters after the last effect
    
%     if preheat_every_effect
        Tf_eph(1:n-1) = Tb(1:n-1) - Delta_Tf_iph;
%     elseif preheat_every_2effects
%         % Odds
%         Tf_eph(i) = Tf_eph(i+1);
%         % Evens
%         Tf_eph(i) = Tb(i) - Delta_Tf_iph;
%     end
    
    % Specific enthalpies difference between feedwater passing through the preheaters
    Hf_Tf_eph = enthalpy_saturated_liquid_water(Tf_eph);
    Delta_H_eph = Hf_Tf_eph(1:n-1) - Hf_Tf_eph(2:n);
    
    %%
    % Specific enthalpies difference between brine and feedwater in the 1st effect
    Hb_Tb(1) = swenthalpy(Tb(1),Pv(1),Xb(1));
    
    % Calculate mass flow rate of feedwater into the 1st effect
    
    
end


%% Simplified MED code - JUST FOR REFERENCE
while(false)
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
end