function [c_t] = find_heatcapacity(Temp, ice, salt, porosity)

%this function calculates the temperature dependent thermal conductivity of
%ice/mixture of ice and salt with the


c_t = (((ice).*(1925 .* (Temp ./ 250)))+ (salt.*2000) + (porosity.*0));


end