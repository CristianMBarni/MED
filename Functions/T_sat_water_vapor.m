function [T] = T_sat_water_vapor(P)

% P is pressure in kPa

T = (42.6776-3892.7/(ln(P/1000)-9.48654)) - 273.15;
end