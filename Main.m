%% Initialization
clc
close all
clearvars
path(pathdef)

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
pre_PlateHTX_flag = 0; % 0/1- Plate heat exchangers of the feedwater entering the down condenser

% Actual inputs
Ms = 1; % Mass flow of motive steam
Q_Loss = 0.01; % Fraction of thermal losses in each effect
Ts_sat = 100; % Temperature of saturated motive steam
Tv1 = 65; % Vapor temperature in the first effect
Tvn = 40; % Vapor temperature in the last effect
Tf_dc_out = 35; % Temperature of feedwater leaving the down condenser
pre_Tsw_out = 30; % Temperature of the seawater leaving the plate HXT
Tsw = 25; % Temperature of seawater
Tv_Loss = 0.2; % Temperature losses by pressure losses (Heat transfer inside each effect is at Tv-Tv_loss)
Xsw = 0.042; % Salinity of seawater in weight %
Xb1 = Xsw+0.003; % Salinity in weight % to be assumed when starting iteration - Calculate mass flow rate of feedwater for the 1st effect
Xbn_max = 0.072; % Maximum allowed salinity of the last effect

% NEEDS REVIEW - El Dessouky, 4.3.2 - last sentence
GOR_guess = 0.98*n; % GOR is aproximately equal to the number of effects
Md = GOR_guess*Ms; % Equation invented randomly, needs reviewing
D = 0;

%% Problem solving
while abs(Md - sum(D)) > 1e-4
    if sum(D) ~= 0
        Md = sum(D);
    end
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
    Delta_T_htx(2:n) = (Tb(2:n)+Tv(2:n))/2; % CHECK THIS EQUATION
    
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
    %         Tf_eph(1:n-1) = Tb(1:n-1) - Delta_Tf_iph;
    %     elseif preheat_every_2effects
    %         % Odds
    %         Tf_eph(i) = Tf_eph(i+1);
    %         % Evens
    %         Tf_eph(i) = Tb(i) - Delta_Tf_iph;
    %     elseif no_preheaters
    Tf_eph = Tf_dc_out;
    %     end
    
    % Specific enthalpies difference between feedwater passing through the preheaters
    %     Hf_Tf_eph = enthalpy_saturated_liquid_water(Tf_eph);
    %     Delta_H_eph = Hf_Tf_eph(1:n-1) - Hf_Tf_eph(2:n);
    
    %%
    % Specific enthalpies difference between brine and feedwater in the 1st effect
    Hb_Tb(1) = swenthalpy(Tb(1),Pv(1),Xb(1));
    Hf_Tf(1) = swenthalpy(pre_Tsw_out,Pv(1),Xb(1));
    Delta_H_iph(1) = Hb_Tb(1) - Hf_Tf(1);
    
    % NEEDS REVIEW
    % Calculate mass flow rate of feedwater into the 1st effect
    if sum(D) == 0
        F_total = Md*Xbn_max/(Xbn_max - Xsw); % Equation invented randomly - CHECK THIS EQUATION
    else
        F_total = D(1)*Xb(1)/(Xb(1) - Xsw);
    end
    %     if parallel_feed
    F(1:n) = F_total/n;
    %     elseif forward_feed
    %
    %     end
    
    %     if flashbox_present
    % Flashing of distillate crossflowing across effects in the distillate boxes
    %
    %     end
    
    %% Calculate effects
    for i = 1:n
        % Calculate Xb(i)
        LHv_evap(i) = latent_heat_water_evaporation(Tv(i));
        Qd_flash(i) = 0;
        if i ~= 1
            Q(i) = (Qv_remain_out(i-1) + Qd_flash(i))*(1-Q_Loss);
            Xb(i) = Xb(i-1);
            
            Hb_Tb(i) = swenthalpy(Tb(i),Pv(i),Xb(i));
            Hf_Tf(i) = swenthalpy(pre_Tsw_out,Pv(i),Xb(i));
            Delta_H_iph(i) = Hb_Tb(i) - Hf_Tf(i);
        end
        
        Q_iter = F(i)*(Delta_H_iph(i) + LHv_evap(i)*(Xb(i)-Xsw)/Xb(i)); % Heat load = (heating feed + evaporating)
        while Q(i) > Q_iter
            Xb(i) = Xb(i) + Xb1*0.00001;
            if Xb(i) > Xbn_max
                error(['Brine salinity got bigger than the maximum allowed on effect ' num2str(i)])
            end
            Hb_Tb(i) = swenthalpy(Tb(i),Pv(i),Xb(i));
            Hf_Tf(i) = swenthalpy(pre_Tsw_out,Pv(i),Xb(i));
            Delta_H_iph(i) = Hb_Tb(i) - Hf_Tf(i);
            Q_iter = F(i)*(Delta_H_iph(i) + LHv_evap(i)*(Xb(i)-Xsw)/Xb(i)); % Heat load = (heating feed + evaporating)
        end
        
        while Q(i) < Q_iter
            Xb(i) = Xb(i) - Xb1*0.00001;
            Hb_Tb(i) = swenthalpy(Tb(i),Pv(i),Xb(i));
            Hf_Tf(i) = swenthalpy(pre_Tsw_out,Pv(i),Xb(i));
            Delta_H_iph(i) = Hb_Tb(i) - Hf_Tf(i);
            if Xb(i) <= Xsw
                error(['Something is wrong with your problem formulation. Your brine salinity in the ' num2str(i) ' effect is lower than the feedwater salinity'])
            end
            Q_iter = F(i)*(Delta_H_iph(i) + LHv_evap(i)*(Xb(i)-Xsw)/Xb(i)); % Heat load = (heating feed + evaporating)
        end
        
        %         if feedwater_preheater
        %
        %         else
        Q_eph(i) = 0;
        %         end
        V_evap(i) = F(i)*((Xb(i)-Xsw)/Xb(i));
        B_evap(i) = F(i) - V_evap(i);
        
        S(i) = (B_evap(i)*Xb(i)*10)/1000;
        if i == 1
            S_out(i) = S(i);
            Xb_out = Xb(i);
        else
            S_out(i) = S(i) + S_out(i-1);
            Xb_out(i) = ((S_out(i-1) + S(i))*1000)/((B_evap(i) + B(i-1))*10);
        end
        
        %         if crossflow
        %
        %         else
        V_b_flash(i) = 0;
        B_b_flash(i) = 0;
        B_b_flash_remain(i) = 0;
        Hb_b_flash_remain(i) = 0;
        Qv_b_flash(i) = 0;
        Tb_b_flash(i) = Tb(i);
        %         end
        
        V(i) = V_evap(i) + V_b_flash(i);
        B(i) = B_evap(i) + B_b_flash(i);
        
        % NO IDEA WHERE THIS EQUATION CAME FROM
        Tb_out(i) = (B_evap(i)*Hb_Tb(i) + B_b_flash_remain(i)*Hb_b_flash_remain(i))/((B_evap(i)*Hb_Tb(i)/Tb(i)) + (B_b_flash_remain(i)*Hb_b_flash_remain(i)/Tb_b_flash(i)));
        
        %         LHv_Tv_out(i) = latent_heat_water_evaporation(Tv_out(i));
        %         V_eph(i) = Q_eph(i)/LHv_Tv_out(i);
        %         V_evap_remain(i) = V_evap(i) - V_eph(i);
        V_evap_remain(i) = V_evap(i);
        
        % p.201
        Qv_evap(i) = V_evap(i)*LHv_evap(i);
        Qv_evap_remain(i) = V_evap_remain(i)*LHv_evap(i);
        
        % Crossflow of distillate mass balance
        if i == 1
            D_created(i) = Ms;
            D(i) = D_created(i);
            D_enter_next(i) = 0;
        else
            D_created(i) = V(i-1);
            D(i) = V(i-1) + D_enter_next(i-1);
            D_enter_next(i) = 0;
        end
        
        % Distillate temperature outlet from effect(i)
        Hv_Tv(i) = saturated_water_enthalpy(Tv(i));
        if i == 1
            Td_out(i) = Ts_sub;
        elseif i == 2
            Td_out(i) = V(i-1)*Hv_Tv(i)/(V(i-1)*Hv_Tv(i)/Tv(i));
            % Td_out(i) is the balance from the below different fluids, something along the lines of the below equation
            % V(i-1); D_d_flash_remain(i); V_d_flash(i); D_d_flash_remain_Ej_s; V_d_flash_Ej_j
            % Tb_out(i) = (B_evap(i)*Hb_Tb(i) + B_b_flash_remain(i)*Hb_b_flash_remain(i))/((B_evap(i)*Hb_Tb(i)/Tb(i)) + (B_b_flash_remain(i)*Hb_b_flash_remain(i)/Tb_b_flash(i)));
        else
            Td_out(i) = V(i-1)*Hv_Tv(i)/(V(i-1)*Hv_Tv(i)/Tv(i));
            % Same as i == 2, but less fluids
            % V(i-1); D_d_flash_remain(i); V_d_flash(i)
        end
        
        Qv_remain_out(i) = Qv_evap_remain(i);
        %         Qd_flash(i) = 0;
        Tv_out(i) = Td_out(i) - Tv_Loss; %INVENTED EQUATION - certainly wrong, needs thinking
    end
    %p.202
    % Rank_Md_return = Ms + E_Ms_Tot;
    % Rank_Td_return = f(D(1);Ts_sub;Td_out(2));
    Rank_Md_return = Ms;
    Rank_Td_return = Ts_sub;
    
    D_Total = sum(D) + V(n) - Rank_Md_return;
    
    % Condenser
    Qd_flash_cond = 0;
    Qv_remain_out(n) = Qv_remain_out(n); % - a bunch of things
    Qdc_Vapor = Qv_remain_out(n) + Qd_flash_cond;
    
    % Preheaters
    Hf_Tf_dc_out = swenthalpy(Tf_dc_out,Pv(1),Xsw);
    Hsw_Tsw = swenthalpy(Tsw,Pv(1),Xsw);
    if Hf_Tf_dc_out - Hsw_Tsw == 0
        Mcw_dc_in = 0;
        Mcw_dc_out_reject = 0;
    else
        Mcw_dc_in = Qdc_Vapor/(Hf_Tf_dc_out - Hsw_Tsw);
        Mcw_dc_out_reject = Mcw_dc_in - F_total;
    end
    Delta_H_dc = Hf_Tf_dc_out - Hsw_Tsw;
    Qdc_Feed = F_total*Delta_H_dc;
    
    Delta_Qdc = Qdc_Vapor - Qdc_Feed;
    Ratio_Qdc = Qdc_Vapor/Qdc_Feed;
    
    % Check condenser operational boundaries
    if Ratio_Qdc < 1.08
        error(['Condensator ratio of heat is too low (' num2str(Ratio_Qdc) ')'])
    end
    if Ratio_Qdc > 2.5
        warning(['Condensator ratio of heat is too high (' num2str(Ratio_Qdc) ')'])
    end
    
    PR = D_Total/Ms;
    
    % Heat transfer coefficients and areas
    U(n) = 1e-3*1939.4 + 1.40562*Tb(n) - 0.0207525*Tb(n)^2 + 0.0023186*Tb(n)^3;
    Ae = Q(n)/(U(n)*(Tv_out(n-1) - Tb(n)));
    U(1) = Q(1)/(Ae*(Ts_sub - Tb(1)));
    
    for i = 2:n-1
        U(i) = Q(i)/(Ae*(Tv_out(i-1) - Tb(i)));
    end
    
    Udc = 1e-3*1617.5 + 0.1537*Tv_out(n) + 0.1825*Tv_out(n)^2 - 0.00008026*Tv_out(n)^3;
    LMTD_dc = (Tf_dc_out - Tsw)/log((Tv_out(n) - Tsw)/(Tv_out(n) - Tf_dc_out));
    Adc = Qdc_Vapor/(Udc*LMTD_dc);
    disp(['Md_guess = ' num2str(Md) ' and sum(D) = ' num2str(sum(D))])
end

GOR = sum(D)/Ms % Gain Output Ratio - also called Performance Ratio (PR)
sA = (Ae*n + Adc)/sum(D) % Specific Heat Transfer Area
RR = sum(D)/sum(F) 

%% Functions used in this program

function [Tbpe] = bpe(T,X)
% T is temperature in ºC
% X is salt weight %
A = 8.325e-2 +1.883e-4.*T +4.02e-6.*T.^2;
B = -7.625e-4 +9.02e-5.*T -5.2e-7.*T.^2;
C = 1.522e-4 -3e-6.*T -3e-8.*T.^2;
Tbpe = A.*X + B.*X.^2+C.*X.^3;
end

function [P] = P_sat_water_vapor(T)

% P is pressure in kPa
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
    aux = aux + f(i).*(0.01.*(T +273.15 -338.15)).^(i-1);
end

P = Pc*exp((Tc./(T+273.15)-1).*aux);
end

function hw = saturated_water_enthalpy(T)
hw = 0.117603379 + 4.20689825*T - 0.000800678287*T^2 + 0.00000895478984*T^3 - 2.61832765E-08*T^4;           %[kJ/kg]
end

function h = swenthalpy(T,P,xs)

Po = 100000;
ws = xs/1000; %kg/kg

% Fitting polynomials through EES for Temperatures between 5-95ºC
% hw = 0.117603379 + 4.20689825*T - 0.000800678287*T^2 + 0.00000895478984*T^3 - 2.61832765E-08*T^4;           %[kJ/kg]
hw = saturated_water_enthalpy(T);
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

function [h] = latent_heat_water_evaporation(T)

% T is temperature in ºC
% range between 0.01-200ºc with 0.017% error from steam tables

% h = 2501.897149 -2.407064037*T +1.192217e-3*T^2 -1.5863e-5*T^3; % Appendix A
h = 2499.5698 - 2.204864*T - 2.304e-3*T^2; % Ch 4.2.2 - Example 1
end