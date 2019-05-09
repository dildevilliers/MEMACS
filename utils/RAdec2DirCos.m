function [u,v,w] = RAdec2DirCos(RA,dec)
% function [u,v,w] = PhTh2DirCos(ph,th)
% Try to maintain the ph angles in the pole
dec(dec == 0 & RA ~= 0) = eps;
dec = pi/2 - dec;
RA = RA;
u = real(sin(dec).*cos(RA));
v = real(sin(dec).*sin(RA));
w = real(cos(dec));
end