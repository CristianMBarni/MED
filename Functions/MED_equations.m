function [D,B,X,T,Tv,hv_vap,A] = MED_equations(n,Md,Mf,Xf,U,deltaT,deltaT_loss,Ts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Tn = T(n);
% Bn = B(n);
% Xn = X(n);

T(1) = Ts - deltaT(1);
for i = 2:n
    T(i) = T(i-1) - deltaT(i);
end

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
for i = 2:n
    B(i) = B(i-1) - D(i);
end

%% Salt concentration profile
X(1) = Xf*Mf/B(1);
for i = 2:n
    X(i) = X(i-1)*B(i-1)/B(i);
end

A(1) = D(1)*hv_vap(1)/(U(1)*(Ts-T(1)));
for i = 2:n
    A(i) = D(i)*hv_vap(i)/(U(i)*(deltaT(i)-deltaT_loss));
end

% check_consistence(Tn,T(n),'Temperature')
% check_consistence(Bn,B(n),'Brine flow rate')
% check_consistence(Xn,X(n),'Brine salinity')
end