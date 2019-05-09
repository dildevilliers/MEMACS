function [u,v,w] = AzAlt2DirCos(az,alt)
% function [u,v,w] = PhTh2DirCos(ph,th)
% Try to maintain the ph angles in the pole

alt = wrap2pi(pi/2 - alt);
az = wrap2pi(-az);
alt(alt == 0 & az ~= 0) = eps;
u = real(sin(alt).*cos(az));
v = real(sin(alt).*sin(az));
w = real(cos(alt));
end