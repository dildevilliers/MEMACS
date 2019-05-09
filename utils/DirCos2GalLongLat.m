function [long,lat] = DirCos2GalLongLat(u,v,w)
% function [long,lat] = DirGalLongLat(u,v,w)
long = real(atan2(v,u));
lat = pi/2 - real(acos(w));
end