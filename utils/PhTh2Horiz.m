function [az,alt] = PhTh2Horiz(ph,th)
% function [az,alt] = PhTh2Horiz(ph,th)

az = wrap2pi(-ph);
alt = wrap2pi(pi/2 - th);
end