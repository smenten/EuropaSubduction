function [k_t] = find_thermalconductivity(Temp, ice, salt, porosity)

%this function calculates the temperature dependent thermal conductivity of
%ice/mixture of ice and salt with the


k_t = (((ice).*(651 ./ Temp))+ (salt.*0.65) + (porosity.*0));

end