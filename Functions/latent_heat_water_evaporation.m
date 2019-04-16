function [h] = latent_heat_water_evaporation(T)

% T is temperature in ºC
% range between 0.01-200ºc with 0.017% error from steam tables

% h = 2501.897149 -2.407064037*T +1.192217e-3*T^2 -1.5863e-5*T^3; % Appendix A
h = 2499.5698 - 2.204864*T - 2.304e-3*T^2; % Ch 4.2.2 - Example 1
end

