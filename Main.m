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
d_Cross_flow_flag = 0; % 0/1 - Cross flow (Parallel)
eph_flag = 0; % 0-No preheaters; 1-PH/NCG; 2-alternating PH/NCG; 3-PH; 4-alternating PH
n_ph_NCG = 0; % Effect which NCG will start to act
pre_PlateHTX_flag = 0; % 0/1- Plate heat exchangers of the feedwater entering the down condenser

% Actual inputs
Ms = 1; % Mass flow of motive steam
Q_Loss = 0.01; % Fraction of thermal losses in each effect
Ts_sat = 70; % Temperature of saturated motive steam
Tv1 = 65; % Vapor temperature in the first effect
Tvn = 40; % Vapor temperature in the last effect
Tf_dc_out = 35; % Temperature of feedwater leaving the down condenser
pre_Tsw_out = 0; % Temperature of the seawater leaving the plate HXT
Tsw = 25; % Temperature of seawater
Tv_Loss = 0.2; % Temperature losses by pressure losses (Heat transfer inside each effect is at Tv-Tv_loss)
% Xsw = 0.042; % Salinity of seawater in weight %
% Xb1 = Xsw+0.003; % Salinity in weight % to be assumed when starting iteration - Calculate mass flow rate of feedwater for the 1st effect
% Xbn_max = 0.072; % Maximum allowed salinity of the last effect

Xsw = 42000; % Salinity of seawater in mg/kg
Xb1 = Xsw + 3000; % Salinity in mg/kg to be assumed when starting iteration - Calculate mass flow rate of feedwater for the 1st effect
Xbn_max = 72000; % Maximum allowed salinity of the last effect

% Inputs added by me
Ts_super = Ts_sat; % Actual temperature of the motive steam - can be bigger in case suerheated steam is used, needs correlation for enthalpy of superheated steam

% NEEDS REVIEW - El Dessouky, 4.3.2 - last sentence
PR_guess = 0.98*n; % GOR is aproximately equal to the number of effects
Md = PR_guess*Ms; % Equation invented randomly, needs reviewing
% Md = 380;
D_Total = 0;
tic
%% Problem solving
while abs(Md - D_Total) > 1e-4
    %     tic
    if D_Total ~= 0
        Md = D_Total;
    end
    % Set/Reset 1st effect salinity
    if not(exist('last_Xb1','var'))
        Xb(1) = Xb1;
    end
    
    %% Vapor and brine temperature profile across effects
    % Set temperature and pressure profile
    Tv(1) = Tv1;
    Tv(n) = Tvn;
    delta_Tv = (Tv(1) - Tv(n))/(n-1);
    Tv = Tv(1):-delta_Tv:Tv(n);
    
    Pv = P_sat_water_vapor(Tv);
    
    % Set temperature profile of the brine
    Tbpe = bpe(Tv,Xb(1)); % Salinity is assumed to be Xb(1) for all effects in order to estimate Tbpe. Error is around 0.12ºC
    Tb = Tv + Tbpe;
    
    % Calculate subcooled temperature of distillate collected from 1st effect
%     Delta_T_htx(2:n) = Tv(1:n-1) - Tb(2:n); % ASSUMED EQUATION, MIGHT BE WRONG
%     Delta_T_htx(1) = sum(Delta_T_htx(2:n))/(n-1); % Estimate that the HTX in the 1st effect happens at the average of this value for the other effects
%     Ts_sub = Tb(1) + Delta_T_htx(1); % Ts_sub is the temperature that the motive steam leaves the 1st effect - Considers subcooled temperature, needs correlation for enthalpy of subcooled water
    
    % MY CHANGE
    Ts_sub = Ts_sat; % For use of only the latent heat
    
    %% If a TVC is present, call TVC subroutines
    
    %% Difference between Tb leaving and Tf entering the effects
    Delta_Tf_iph = Tb(n) - Tf_dc_out;
    
    %% Temperature profile of feedwater input across effects
    Tf_eph(n) = Tf_dc_out; % There are no preheaters between down condenser and last effect
    
    if eph_flag == 0 % No preheaters
        Tf_eph(1:n-1) = Tf_dc_out;
    elseif ismember(eph_flag,[1 3]) % Preheaters every effect
        Tf_eph(1:n-1) = Tb(1:n-1) - Delta_Tf_iph;
        
        if eph_flag == 1 % NCG ejectors are present
            
        end
        
    elseif ismember(eph_flag,[2 4]) % Preheaters every two effects
        Tf_eph(logical(mod(1:n-1,2))) = Tf_eph(i+1); % Odd effects have preheaters
        Tf_eph(not(mod(1:n-1,2))) = Tb(i) - Delta_Tf_iph; % Even effects doesn't have preheaters
        
        if eph_flag == 2 % NCG ejectors are present
            
        end
    end
    
    % Specific enthalpies difference between feedwater passing through the preheaters
    Hf_Tf_eph = seawater_enthalpy(Tf_eph,Pv(1),Xsw);
    Delta_H_eph = Hf_Tf_eph(1:n-1) - Hf_Tf_eph(2:n);
    
    %% Specific enthalpies difference between brine and feedwater in the 1st effect - not necessaty for PF
    %     Hb_Tb(1) = seawater_enthalpy(Tb(1),Pv(1),Xb(1));
    Hf_Tf = Hf_Tf_eph;
    %     Delta_H_iph(1) = Hb_Tb(1) - Hf_Tf(1);
    
    %% Calculate mass flow rate of feedwater into the 1st effect
    % NO EQUATIONS AVAILABLE on Casimiro Thesis
    %     Mf = Md + Mb;
    %     Mb = Mf*Xsw/Xbn_max;
    %     Mf = Md + Mf*Xsw/Xbn_max;
    %     Mf*(1-Xsw/Xbn_max) = Md
    %     Mf = Md/(1-Xsw/Xbn_max);
    
    if CT == 7 || CT == 8 % Low Temperature Parallel Feed || TVC Parallel Feed
        Mf = Md*Xbn_max/(Xbn_max - Xsw/0.7); % Estimated total feedwater mass flow rate for reaching Xbn_max
        F(1:n) = Mf/n; % Feedwater mass flow rate considered equal for each effect
        Xf(1:n) = Xsw;
    elseif CT == 9 || CT == 10 % Low Temperature Forward Feed || TVC Forward Feed
        
    end
    
    %% If NCG steam ejectors are present, call NCG subroutine
    
    %% Calculate 1:n effects
    %     toc
    %     tic
    for i = 1:n
        %% If there is cross flow of distillate between effects, calculates flashbox(i)
        if d_Cross_flow_flag
            
        else
            Qd_flash(i) = 0;
        end
        
        %% Power available for effect(i)
        if i == 1
            %             Hs_super = superheated_water_enthalpy(Ts);
            %             Hs_sub = subcooled_water_enthalpy(Ts_sub);
            %             Q(1) = Ms*(Hs_super - Hs_sub);
            
            % MY CHANGE
            delta_Hs = latent_heat_water_evaporation(Ts_sat);
            Q(1) = Ms*(delta_Hs);
        else
            Q(i) = (Qv_remain_out(i-1) + Qd_flash(i))*(1-Q_Loss);
        end
        
        if CT == 7 || CT == 8 % Low Temperature Parallel Feed || TVC Parallel Feed
            %% Calculate Xb(i)
            if i ~= 1
                Xb(i) = Xb(i-1); % Initial guess
            end
            
            Hb_Tb(i) = seawater_enthalpy(Tb(i),Pv(i),Xb(i));
            Delta_H_iph(i) = Hb_Tb(i) - Hf_Tf(i);
            LHv_evap(i) = latent_heat_water_evaporation(Tv(i));
            
            D_evap(i) = F(i)*(Xb(i)-Xf(i))/Xb(i);
            Q_iter = F(i)*Delta_H_iph(i) + D_evap(i)*LHv_evap(i); % Heat load = (heating feed + evaporating)
            
            while Q(i) > Q_iter % Increases Xb(i) until the power available is equal to the power necessary
                Xb(i) = Xb(i) + Xb1*1e-5;
                if Xb(i) > Xbn_max
                    warning(['Brine salinity got bigger than the maximum allowed on effect ' num2str(i)])
                end
                Hb_Tb(i) = seawater_enthalpy(Tb(i),Pv(i),Xb(i));
                Hf_Tf(i) = seawater_enthalpy(pre_Tsw_out,Pv(i),Xb(i));
                Delta_H_iph(i) = Hb_Tb(i) - Hf_Tf(i);
                Q_iter = F(i)*(Delta_H_iph(i) + LHv_evap(i)*(Xb(i)-Xsw)/Xb(i)); % Heat load = (heating feed + evaporating)
            end
            
            while Q(i) < Q_iter
                Xb(i) = Xb(i) - Xb1*1e-5;
                Hb_Tb(i) = seawater_enthalpy(Tb(i),Pv(i),Xb(i));
                Hf_Tf(i) = seawater_enthalpy(pre_Tsw_out,Pv(i),Xb(i));
                Delta_H_iph(i) = Hb_Tb(i) - Hf_Tf(i);
                if Xb(i) <= Xsw
                    error(['Something is wrong with your problem formulation. Your brine salinity in the ' num2str(i) ' effect is lower than the feedwater salinity'])
                end
                Q_iter = F(i)*(Delta_H_iph(i) + LHv_evap(i)*(Xb(i)-Xsw)/Xb(i)); % Heat load = (heating feed + evaporating)
            end
            
            %% Thermal load used in the external feedwater preheaters
            % If preheater is present in effect(i)
            if (ismember(eph_flag,[1 3]) || (ismember(eph_flag,[2 4]) && mod(i,2))) && not(i == n)
                Q_eph(i) = sum(F(1:i))*Delta_H_eph(i);
            else
                Q_eph(i) = 0;
            end
            
            %% Mass flow rate of vapor, brine and salt produced by evaporation in effect(i)
            V_evap(i) = F(i)*((Xb(i)-Xf(i))/Xb(i));
            
            B_evap(i) = F(i) - V_evap(i);
            
            S(i) = B_evap(i)*Xb(i); % Mass of salt is given in mg
            
            %% Total mass of salt and salinity exiting effect (i)
            if i == 1
                S_out(1) = S(1);
                Xb_out(1) = Xb(1);
            else
                S_out(i) = S(i) + S_out(i-1);
                Xb_out(i) = S_out(i)/(B_evap(i) + B(i-1));
            end
            
            %% Flashing of brine crossflowing across effects (i-1) and (i)
            if i == 1 % No flashing of brine inside the first effect
                V_b_flash(i) = 0;
                B_b_flash_remain(i) = 0;
                Hb_b_flash_remain(i) = 0;
                Tb_b_flash = Tb(i);
                Qv_b_flash(i) = 0;
            else
                % Non-Equilibrium allowance between the hotter and colder brines
                Delta_T_b_NEA(i) = nea(Tb_out(i-1)-Tb(i),Tv(i));
                
                % Final temperature of the brine after flashing
                Tb_b_flash(i) = Tb(i) + Delta_T_b_NEA(i);
                
                % Temperature of vapour produced by flashing
                Tv_b_flash(i) = Tb_b_flash(i) - bpe(Tb_b_flash(i),Xb_out(i));
                
                % Mass flow of vapour produced inside the effect by flash of the brine
                Hb_Tb_out(i-1) = seawater_enthalpy(Tb_out(i-1),Pv(i-1),Xb(i-1));
                Hb_b_flash_remain(i) = seawater_enthalpy(Tb_b_flash(i),Pv(i),Xb(i));
                LHv_b_flash(i) = latent_heat_water_evaporation(Tv_b_flash(i));
                V_b_flash(i) = B(i-1)*(Hb_Tb_out(i-1) - Hb_b_flash_remain(i))/LHv_b_flash(i);
                
                % Mass flow of brine that remains after flash
                B_b_flash_remain(i) = B(i-1) - V_b_flash(i);
                
                % Thermal load released from flashing
                Qv_b_flash(i) = V_b_flash(i)*LHv_b_flash(i);
            end
            
            %% Total mass flow rate of vapor and brine formed inside effect(i)
            V(i) = V_evap(i) + V_b_flash(i);
            B(i) = B_evap(i) + B_b_flash_remain(i);
            
            %% Outlet brine temperature
            % NO IDEA WHERE THIS EQUATION CAME FROM, seems like a weighted mean
            Tb_out(i) = (B_evap(i)*Hb_Tb(i) + B_b_flash_remain(i)*Hb_b_flash_remain(i))/((B_evap(i)*Hb_Tb(i)/Tb(i)) + (B_b_flash_remain(i)*Hb_b_flash_remain(i)/Tb_b_flash(i)));
            
        elseif CT == 9 || CT == 10 % Low Temperature Forward Feed || TVC Forward Feed
            
        end
        
        %% Separation of mass flow rate of vapor formed inside effect(i)
        Tv_out(i) = Tv(i) - Tv_Loss;
        LHv_Tv_out(i) = latent_heat_water_evaporation(Tv_out(i));
        V_eph(i) = Q_eph(i)/LHv_Tv_out(i);
        V_evap_remain(i) = V_evap(i) - V_eph(i);
        
        %% Heat load of vapor formed by evaporation inside effect(i) and amount used to power the next effect
        Qv_evap(i) = V_evap(i)*LHv_evap(i);
        Qv_evap_remain(i) = V_evap_remain(i)*LHv_evap(i);
        
        %% Crossflow of distillate mass balance
        if i == 1
            D_created(1) = Ms; % Distillate produced is equal to the motive steam entering the MED
            D(1) = D_created(1);
            if ismember(CT,[8 10]) % If TVC is present
                if d_Cross_flow_flag
                    
                else
                    D_enter_next(1) = 0;
                end
            else
                D_enter_next(1) = 0;
            end
        else
            D_created(i) = V(i-1);
            D(i) = V(i-1) + D_enter_next(i-1);
            
            % If NCG steam ejectors are present some code needs to be added
            
            % Mass flow of distillate entering next effect
            if d_Cross_flow_flag
                D_enter_next(i) = D(i);
            else
                D_enter_next(i) = 0;
            end
            
        end
        
        %% Distillate temperature outlet from effect(i)
        if i == 1
            Td_out(1) = Ts_sub;
        else
            if ismember(CT,[8 10]) % If TVC is present
                
            else
                % NO IDEA WHERE THIS EQUATIONS CAME FROM, seems like a weighted mean
                if i == 2
                    Td_out(i) = V(i-1)*Hd_Tv_out(i-1)/(V(i-1)*Hd_Tv_out(i-1)/Tv_out(i-1));
                    
                    % V(i-1) is the Vapor formed inside effect (i-1);
                    % D_d_flash_remain(i) is the distillate that remains after the crossflowing distillate suffers flashing;
                    % V_d_flash(i) is the vapor that is produced after the crossflowing distillate suffers flashing;
                    % D_d_flash_remain_Ej_s & V_d_flash_Ej_j are related to the NCG steam ejectors
                else
                    Td_out(i) = V(i-1)*Hd_Tv_out(i-1)/(V(i-1)*Hd_Tv_out(i-1)/Tv_out(i-1));
                    % Same as i == 2, but no NCG fluids
                    % Td_out(i) = (V(i-1)*Hd_Tv_out(i-1) + D_d_flash_remain(i)*Hd_d_flash_remain(i) + V_d_flash(i)*Hd_Tv_out(i-1))/((V(i-1)*Hd_Tv_out(i-1)/Tv_out(i-1)) + (D_d_flash_remain(i)*Hd_d_flash_remain(i)/Tv_d_flash(i)) + (V_d_flash(i)*Hd_Tv_out(i-1)/Tv_out(i-1)));
                    
                end
            end
        end
        
        % Calculating necessary enthalpies
        Hd_Tv_out(i) = saturated_liquid_water_enthalpy(Tv_out(i));
        Hd_Td_out(i) = saturated_liquid_water_enthalpy(Td_out(i));
        
        if ismember(CT,[7 8]) % If Parallel Feed
            Qv_remain_out(i) = (V_evap_remain(i) + V_b_flash(i))*LHv_Tv_out(i);
            Qv(i) = Qv_evap(i) + Qv_b_flash(i);
        else
            Qv_remain_out(i) = (V_evap_remain(i) + V_b_flash(i) + V_f_flash(i))*LHv_Tv_out(i);
            Qv(i) = Qv_evap(i) + Qv_b_flash(i) + Qv_f_flash(i);
        end
    end
    %     toc
    % p.202
    %     tic
    %% Mass flow of condensate returning to the Rankine cycle or Boiler
    % If NCG is not present:
    E_Mm_Tot = 0;
    
    Rank_Md_return = Ms + E_Mm_Tot;
    
    if ismember(CT,[8 10]) % If TVC is present
        % If NCG is present:
        % Rank_Td_return = E_Tcw_out_2;
        % Else
        % Rank_Td_return = Ts_sub;
        % End
    else
        Hs_Ts_sub = saturated_liquid_water_enthalpy(Ts_sub); % If Ts_sub is subcooled, change correlation
        Rank_Td_return = (D(1)*Hd_Tv_out(1) + Ms*Hs_Ts_sub + D(2)*Hd_Td_out(2))/(D(1)*Hd_Tv_out(1)/Tv_out(1) + Ms*Hs_Ts_sub/Ts_sub + 1/Td_out(2));
    end
    
    %% Total mass flow of Distillate produces
    if d_Cross_flow_flag
        D_Total = D(n) + V(n);
    else
        D_Total = sum(D) + V(n) - Rank_Md_return;
    end
    
    %% Flashing of Distillate inside down condenser
    % NEA between hotter distillate entering the dc and colder after flashing
    Delta_T_d_NEA_cond = nea(Td_out(n)-Tv_out(n),Tv_out(n));
    
    % Temperature of vapor produced with the flashing
    Tv_d_flash_cond = Tv_out(n) + Delta_T_d_NEA_cond;
    
    % Mass flow of vapor produced through flash
    Hd_Td_out(n) = saturated_liquid_water_enthalpy(Td_out(n));
    Hd_d_flash_remain_cond = saturated_liquid_water_enthalpy(Tv_d_flash_cond);
    LHv_d_flash_cond = latent_heat_water_evaporation(Tv_d_flash_cond);
    V_d_flash_cond = D(n)*(Hd_Td_out(n) - Hd_d_flash_remain_cond)/LHv_d_flash_cond;
    
    % Mass flow of distillate remaining
    D_d_flash_remain_cond = D(n) - V_d_flash_cond;
    
    % Thermal load released from flashing
    Qd_flash_cond = V_d_flash_cond*(LHv_d_flash_cond - Hd_Tv_out(n));
    
    %% Heat load actually entering the down condenser
    if ismember(CT,[8 10]) && false % If TVC is present && some conditions
        error('Wrong')
    else
        if not(ismember(CT,[8 10])) % If TVC is NOT present
            TVC_Me_1 = 0;
        end
        
        Hv_Tv = 1; % Only used when TVC is present, needs to be adressed before using
        Hv_b_flash = 1; % Only used when TVC is present, needs to be adressed before using
        if CT == 8 % TVC Parallel Feed
            Qv_remain_out(n) = Qv_remain_out(n) - (Hv_Tv(n)*TVC_Me_1*V_evap(n)/(V_evap(n) + V_b_flash(n))) - (Hv_b_flash(n)*TVC_Me_1*V_b_flash(n)/(V_evap(n) + V_b_flash(n)));
        elseif CT == 10 % TVC Forward Feed
            Qv_remain_out(n) = Qv_remain_out(n) - Hv_Tv(n)*TVC_Me_1;
        end
        
        E_Me_Vap_1 = 0; % If NCG is not present
        Qdc_Vapor = Qv_remain_out(n) + Qd_flash_cond - E_Me_Vap_1;
    end
    
    %p.203
    %% Flat Plate cooling water preheaters
    if pre_PlateHTX_flag % If Plate Heat Exchangers installed
        
    else
        Hf_Tf_dc_out = seawater_enthalpy(Tf_dc_out,Pv(1),Xsw);
        Hsw_Tsw = seawater_enthalpy(Tsw,Pv(1),Xsw);
        if Hf_Tf_dc_out - Hsw_Tsw == 0
            Mcw_dc_in = 0;
            Mcw_dc_out_reject = 0;
        else
            Mcw_dc_in = Qdc_Vapor/(Hf_Tf_dc_out - Hsw_Tsw);
            Mcw_dc_out_reject = Mcw_dc_in - Mf;
        end
        Delta_H_dc = Hf_Tf_dc_out - Hsw_Tsw;
        Qdc_Feed = Mf*Delta_H_dc;
    end
    
    % Difference and ratio between thermal load passed into the feedwater
    % and released by the vapor entering and formed at the down condenser
    Delta_Qdc = Qdc_Vapor - Qdc_Feed;
    Ratio_Qdc = Qdc_Vapor/Qdc_Feed;
    
    %% Check down condenser operational boundaries
    %% If heatload is too low
    if Ratio_Qdc < 1.08
        if CT == 7 || CT == 8 % Low Temperature Parallel Feed || TVC Parallel Feed
            error(['Condensator ratio of heat is too low (' num2str(Ratio_Qdc) ')'])
        elseif CT == 9 || CT == 10 % Low Temperature Forward Feed || TVC Forward Feed
            if Xb_out(n) <= Xbn_max
                if exist('last_Xb1','var') && Xb(1) == last_Xb1
                    warning('Maybe the change in Xb(1) implemented is wrong, check it')
                    error(['Condensator ratio of heat is too low (' num2str(Ratio_Qdc) ') and changing the brine salinity does not solves the problem'])
                elseif Delta_Tf_iph > Tv(1) - Tv(n)
                    error(['Condensator ratio of heat is too low (' num2str(Ratio_Qdc) ') and Delta_Tf_iph > Tv(1) - Tv(n)'])
                end
                last_Xb1 = Xb(1);
                Xb(1) = Xb1*1.01;
                continue
            else
                error(['Condensator ratio of heat is too low (' num2str(Ratio_Qdc) ') and Brine is already over maximum salinity'])
            end
        end
    end
    
    %% If heatload is too high
    if ismember(CT,[8 10])
        if Ratio_Qdc > 1.18
            warning(['Condensator ratio of heat is too high (' num2str(Ratio_Qdc) ')'])
        end
    else
        if Ratio_Qdc > 2.5
            if ismember(CT,[9 10]) && (Xb_out(n) <= Xb_out(n-1) || Xb_out(n) >= Xbn_max || Xb_out(n) <= Xf(1))
                if exist('last_Xb1','var') && Xb(1) == last_Xb1
                    warning('Maybe the change in Xb(1) implemented is wrong, check it')
                    error(['Condensator ratio of heat is too low (' num2str(Ratio_Qdc) ') and changing the brine salinity does not solves the problem'])
                end
                last_Xb1 = Xb(1);
                Xb(1) = Xb1*0.99;
                continue
            else
                warning(['Condensator ratio of heat is too high (' num2str(Ratio_Qdc) ')'])
            end
        end
    end
    
    %% Performance Ratio (PR)-  also called Gain Output Ratio (GOR)
    PR = D_Total/Ms;
    
    %% Heat transfer coefficients and areas of effects
    U(n) = ue(Tb(n));
    Ae = Q(n)/(U(n)*(Tv_out(n-1) - Tb(n)));
    U(1) = Q(1)/(Ae*((Ts_super + Ts_sub)/2 - Tb(1)));
    
    for i = 2:n-1
        U(i) = Q(i)/(Ae*(Tv_out(i-1) - Tb(i)));
    end
    
    %% Heat transfer coefficient and area of down condenser
    Udc = uc(Tv_out(n));
    if pre_PlateHTX_flag
        LMTD_dc = (Tf_dc_out - pre_Tsw_out)/log((Tv_out(n) - pre_Tsw_out)/(Tv_out(n) - Tf_dc_out));
    else
        LMTD_dc = (Tf_dc_out - Tsw)/log((Tv_out(n) - Tsw)/(Tv_out(n) - Tf_dc_out));
    end
    Adc = Qdc_Vapor/(Udc*LMTD_dc);
    %     toc
    disp(['Md_guess = ' num2str(Md) ' and D_Total = ' num2str(D_Total)])
end

PR = D_Total/Ms % Gain Output Ratio - also called Performance Ratio (PR)
sA = (Ae*n + Adc)/D_Total % Specific Heat Transfer Area
RR = D_Total/Mf
toc

%% Pure Water Thermodynamic Properties
function hw = saturated_liquid_water_enthalpy(T)
% Appendix A.5 from El Dessouky
% 5 < T < 200 ºC
% hw is the Enthalpy of saturated liquid water in kJ/kg

hw = -0.033635409 + 4.207557011*T - 6.200339e-4*T.^2 + 4.459374e-6*T.^3; % El Dessouky
end

function hw = saturated_water_vapor_enthalpy(T)
% Appendix A.6 from El Dessouky
% 0.01 < T < 200 ºC
% hw is the enthalpy of saturated water vapor in kJ/kg

hw = 2501.689845 + 1.806916015.*T + 5.087717e-4.*T.^2 -1.1221e-5.*T.^3; % El Dessouky
end

function h = latent_heat_water_evaporation(T)
% Appendix A.7 from El Dessouky
% 0.01 < T < 200 ºC
% h is the latent heat of water evaporation in kJ/kg

h = 2501.897149 -2.407064037*T +1.192217e-3*T^2 -1.5863e-5*T^3; % Appendix A
end

function P = P_sat_water_vapor(T)
% Appendix A.10 from El Dessouky
% P is pressure in kPa
% 5 < T < 200 ºC

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

P = Pc.*exp((Tc./(T+273.15)-1).*aux);
end

%% Seawater Thermodynamic Properties
function h_sw = seawater_enthalpy(T,P,X)
% Eq 25 on Table 9 of Nayar, et. al. - https://doi.org/10.1016/j.desal.2016.02.024
% 10 < T < 120ºC
% 0 < P < 12000 kPa
% 0 < X < 120000 ppm
% h_sw is enthalpy in kJ/kg

% Following commented code improves results a little, but slows down the program
% Po = P_sat_water_vapor(T); % Reference pressure in kPa
% h_w = saturated_liquid_water_enthalpy(T)*1000; % Water enthalpy at pressure P0 in J/kg

Po = 101; % Reference pressure in kPa
h_w = 141.355 + 4202.07*T -0.535*T.*T + 0.004*T.*T.*T; % Water enthalpy at pressure P0 in J/kg
S = X./1e3; % Salinity in g/kg
S_weight = X./1e6; %Salinity in kg/kg

b1 = -2.34825e4;
b2 = 3.15183e5;
b3 = 2.80269e6;
b4 = -1.44606e7;
b5 = 7.82607e3;
b6 = -4.41733e1;
b7 = 2.1394e-1;
b8 = -1.99108e4;
b9 = 2.77846e4;
b10 = 9.72801e1;
h_sw_Po = h_w - S_weight.*(b1 + b2.*S_weight + b3.*S_weight.^2 + b4.*S_weight.^3 + b5.*T + b6.*T.^2 + b7.*T.^3 + b8.*S_weight.*T + b9.*S_weight.^2.*T + b10.*S_weight.*T.^2);
% Seawater entahlpy at pressure Po in J/kg

a1 = 996.7767;
a2 = -3.2406;
a3 = 0.0127;
a4 = -4.7723e-5;
a5 = -1.1748;
a6 = 0.01169;
a7 = -2.6185e-5;
a8 = 7.0661e-8;
h_sw = 1e-3.*(h_sw_Po + 1e-3.*(P-Po).*(a1 + a2.*T + a3.*T.^2 + a4.*T.^3 + S.*(a5 + a6.*T + a7.*T.^2 + a8.*T.^3)));
% Seawater entahlpy at pressure P in kJ/kg
end

function rhow = seawater_density(T,X)
% Appendix A.1 from El Dessouky
% 10 < T < 180 ºC
% 0 < X < 160000 ppm
% rhow is the seawater density in kg/m^3

B = (2*X/1000-150)/150;
G1 = 0.5;
G2 = B;
G3 = 2*B*B - 1;

A1 = 4.032219*G1 + 0.115313*G2 + 3.26e-4*G3;
A2 = -0.10819*G1 + 1.571e-3*G2 - 4.23e-4*G3;
A3 = -0.012247*G1 + 1.74e-3*G2 -9e-6*G3;
A4 = 6.92e-4*G1 - 8.7e-5*G2 - 3.5e-5*G3;

A = (2*T-200)/110;

F1 = 0.5;
F2 = A;
F3 = 2*A*A - 1;
F4 = 4*A*A*A - 3*A;

rhow = 1e3*(A1*F1 + A2*F2 + A3*F3 + A4*F4);
end

function cp = seawater_cp(T,X)
% Appendix A.2 from El Dessouky
% 20 < T < 180 ºC
% 20000 < X < 160000 ppm
% Cp is the seawater specific heat at constant pressure in kJ/kgºC

A = 4206.8 - 6.6197*X + 1.2288e-2*X^2;
B = -1.1262 + 5.4178e-2*X - 2.2719e-4*X^2;
C = 1.2026e-2 - 5.3566e-4*X + 1.8906e-6*X^2;
D = 6.8777e-7 + 1.517e-6*X - 4.4268e-9*X^2;

cp = (A + B*T + C*T^2 + D*T^3)*1E-3;
end

%% Thermodynamic Losses
function Tbpe = bpe(T,X)
% Appendix B.1 from El Dessouky
% 10 < T < 180 ºC
% 10000 < x < 160000 ppm

X_weight = X/1e4; % X in weight percentage, 1 < X_weight < 16 %

A = 8.325e-2 +1.883e-4.*T +4.02e-6.*T.^2;
B = -7.625e-4 +9.02e-5.*T -5.2e-7.*T.^2;
C = 1.522e-4 -3e-6.*T -3e-8.*T.^2;
Tbpe = A.*X_weight + B.*X_weight.^2 + C.*X_weight.^3;
end

function Tnea = nea(deltaTb,Tv)
% Appendix B.2 from El Dessouky
% Tnea is the non-equilibrium allowance in MEE - Miyatake et al. (1973)
% All temperatures are in ºC
% deltaTb is the temperature difference of boiling brine in effects j and j-1
% deltaTb = Tb(j-1) - Tb(j)
% Tv is Tv(j) - Vapor temperature in effect j
% Tv(j) = T(j) - Tbpe(j)

Tnea = 33*deltaTb^0.55/Tv;
end

%% Heat Transfer Coefficients
function u = ue(T)
% Appendix C.6 from El Dessouky
% Overall heat transfer coefficient for a fouled evaporator
% Water forming a falling film on a horizontal tube bundle and vapor condensing inside the tubes
% 2 < u < 4 kW/m^2ºC
% T is the evaporation temperature in ºC

u = 1e-3.*(1939.4 + 1.40562.*T -0.0207525.*T.^2 + 0.0023186.*T.^3);
end

function u = uc(T)
% Appendix C.6 from El Dessouky
% Overall heat transfer coefficient for a fouled condenser
% Vapor condensing on the outside surface and seawater flow on the tube side
% 2 < u < 4 kW/m^2ºC
% T is the condensation temperature in ºC

u = 1e-3.*(1617.5 + 0.1537.*T + 0.1825.*T.^2 -0.00008026.*T.^3);
end