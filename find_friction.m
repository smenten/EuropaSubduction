function [miu] = find_friction (Temp, addsalt)


if Temp <= 213 & addsalt == 0

    miu = 0.7;

end

% if Temp > 213 & addsalt = 0
% 
%     miu = 0.
% 
% end

if Temp <= 213 & addsalt == 1

    miu = 0.7;

end


% if Temp > 213 & addsalt = 1
% 
%     miu  = 
% 
% end

end