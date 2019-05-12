function [long,lat] = RAdec2GalLongLat(RA,dec)
% function [long,lat] = RAdec2GalLongLat(RA,dec)
% Calculates equitorial coordinates RA and dec from galactic coordinates.
% Inputs
% - long: Galactic longitude (radians)
% - lat: Galactic latitude (radians)
% Outputs
% - RA: right ascension (radians)
% - dec: declination (radians)

[galCoords] = celestial.coo.coco([RA dec],'j2000.0','g','r','r');
long = wrap2pi(galCoords(:,1)); 
lat = galCoords(:,2); 
end