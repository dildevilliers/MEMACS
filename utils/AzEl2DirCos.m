function [u,v,w] = AzEl2DirCos(az,el)
% function [u,v,w] = AzEl2DirCos(az,el)
% Try to maintain the az angles in the pole
el(el == pi/2 & az ~= 0) = pi/2+eps;
u = real(sin(az).*cos(el));
v = real(sin(el));
w = real(cos(az).*cos(el));
end