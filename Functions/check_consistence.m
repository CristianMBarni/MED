function [] = check_consistence(val1,val2,prop_name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if abs(val1 - val2) > 0.02
    warning(['Specified final ' prop_name ' is different from calculated by ' num2str(abs(val1 - val2))])
end

end

