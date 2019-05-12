function [RA,dec] = GalLongLat2RAdec(long,lat)
% function [RA,dec] = GalLongLat2RAdec(long,lat)
% Calculates equitorial coordinates RA and dec from galactic coordinates.
% Inputs
% - long: Galactic longitude (radians)
% - lat: Galactic latitude (radians)
% Outputs
% - RA: right ascension (radians)
% - dec: declination (radians)

[equCoords] = celestial.coo.coco([long lat],'g','j2000.0','r','r');
RA = wrap2pi(equCoords(:,1)); %right ascension
dec = wrap2pi(equCoords(:,2)); %declination
end