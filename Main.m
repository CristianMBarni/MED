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
% Parameters:
n = 12;
T(n) = 41;
% deltaTbpe_loss
Tf = 35;
Tol = 1e-4;
% U(1)
X(n) = 72;

% Inputs:
Ms = 12; 
Ts = 74;
Tsw = 25;
Xf = 35;
% Qextra_loss
% deltaT_preheat_loss
U_reduction = 0.95;

% First case
% n = 6; % Number of effects
% Ts = 100; % ºC
% Md = 1; % kg/s
% Xf = 42000; % ppm or mg/kg
% X(n) = 70000; % ppm or mg/kg
% T(n) = 40; % ºC
% deltaT_loss = 2; % ºC
% Tf = 35; % ºC
% Tcw = 25; % ºC

% Second case
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

% Initial guesses
Mf = 50;

Areas = zeros(1,n);
first_iteration = true;
iter = 0;
while first_iteration || (any(abs((A(2:n) - A(1:n-1))) > Tol) || any(abs((A - Areas)) > 1e-2))
    %% Iteration setup
    if first_iteration
        first_iteration = false;
    else
        Areas = A;
    end

    %% Heat transfer coefficients
    U(1) = ue(Ts);
    for i = 2:n
        U(i) = U_reduction*U(i-1);
    end

    %% Initial temperature profile
    deltaT_total = Ts-T(n); % Overall temperature difference
    
    deltaT(1) = deltaT_total/(U(1)*sum(1./U));
    for i = 2:n
        deltaT(i) = deltaT(1)*U(1)/U(i);
    end
    
    T(1) = Ts - deltaT(1);
    for i = 2:n
        T(i) = T(i-1) - deltaT(i);
    end
    
    for i = 1:n
        hv_vap(i) = latent_heat_water_evaporation(T(i));
    end
    
    %% Distillate flow rate
    hs_vap = latent_heat_water_evaporation(Ts);
    D(1) = Ms*hs_vap/hv_vap(1);
    
    Md = 0;
    for i = 1:n
        Md = Md + D(1)*hv_vap(1)/hv_vap(i);
    end
    
    aux = 0;
    for i = 1:n
        aux = aux + hv_vap(1)/hv_vap(i);
    end
    D(1) = Md/aux;
    
    for i = 2:n
        D(i) = D(1)*hv_vap(1)/hv_vap(i);
    end
    
    %% Brine flow rate
    
    B(n) = (Xf/(X(n)-Xf))*Md; % Mass of brine released in the last efect
    
    B(1) = Mf - D(1);
    for i = 2:n
        B(i) = B(i-1) - D(i);
    end
    
    %% Salt concentration profile
    X(1) = Xf*Mf/B(1);
    for i = 2:n
        X(i) = X(i-1)*B(i-1)/B(i);
    end
    
    %% Areas calculation
    for i = 1:n
        deltaTbpe_loss = bpe(T(i),X(i));
        A(i) = D(i)*hv_vap(i)/(U(i)*(deltaT(i)-deltaTbpe_loss));
    end
    iter = iter + 1;
    disp(['Iteration ' num2str(iter)])
    abs((A(2:n) - A(1:n-1)))
end

%%
B(n) = Xf/(X(n)-Xf)*Md;

for i = 2:n
    Q(i) = Q(1);
end

deltaT_total = Ts - T(n);
deltaT(1) = deltaT_total/(U(1)*sum(1./U));
for i = 2:n
    deltaT(i) = deltaT(1)*U(1)/U(i);
    T(i) = T(i-1) - deltaT(i);
end

j = 0;
D1 = 12;
D(1) = 1;
while abs(D1-D(1)) > 1e-3
    D1 = D(1);
    for i = 2:n
        D(i) = D(1)*hv_vap(1)/hv_vap(i);
    end
    Md = 0;
    for i = 1:n
        Md = Md + D(i)*hv_vap(1)/hv_vap(i);
    end
    D(1) = Md/sum(hv_vap(1)./hv_vap(1:n));
    j = j + 1;
end

B(1) = Mf - D(1);
for i = 2:n
    B(i) = B(i-1) - D(i);
end

X(i) = X(i-1)*B(i-1)/B(i);

A(i) = D(i)*hv_vap(i)/(U(i)*(deltaT(i) - deltaTbpe_loss));