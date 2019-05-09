function [u,v,w] = ElAz2DirCos(ep,al)
% function [u,v,w] = ElAz2DirCos(ep,al)
% Try to maintain the ep angles in the pole
al(al == pi/2 & ep ~= 0) = pi/2+eps;
u = real(sin(al));
v = real(cos(al).*sin(ep));
w = real(cos(al).*cos(ep));
end