function [ph,th] = Horiz2PhTh(az,alt)
% function [ph,th] = Horiz2PhTh(az,alt)
% Converts horizontal coordinates (az from north and alt from horison) to 
% local [ph,th] coordinates.
ph = wrap2pi(-az);
th = wrap2pi(pi/2 - alt);
end