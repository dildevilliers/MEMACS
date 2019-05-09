function [x,y] = DirCos2Mollweide(u,v,w)
% function [x,y] = DirCos2Mollweide(u,v,w)

long = -real(atan2(v,u));
lat = pi/2 - real(acos(w));

% R is really R.*S (R-radius, S-scale factor)
R = 1;
x = 2.*long.*sqrt(2).*R.*cos(lat)./pi;
y = sqrt(2).*R.*sin(lat);

end