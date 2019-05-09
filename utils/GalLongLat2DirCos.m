function [u,v,w] = GalLongLat2DirCos(long,lat)
% function [u,v,w] = PhTh2DirCos(ph,th)
% Try to maintain the ph angles in the pole
lat(lat == 0 & long ~= 0) = eps;
lat = pi/2 - lat;
long = long;
u = real(sin(lat).*cos(long));
v = real(sin(lat).*sin(long));
w = real(cos(lat));
end