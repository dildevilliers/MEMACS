function [long,lat] = Horiz2GalLongLat(az,alt,julDate,location)
% function [galLong,galLat] = Horiz2Gal(az,alt,julDate,location)
% Calculates galactic coordinates from horizontal coordinates.
% Inputs
% - az: Horizontal azimuth (radians)
% - alt: Horizontal altitude/elevation (radians)
% - julDate: The Julian date (nr of days)
% - location: earth longitude and latitude in radians [long, lat]
% Outputs
% - long: Galactic longitude (radians)
% - lat: Galactic latitude (radians)

[RA,dec] = Horiz2RAdec(az,alt,julDate,location);
% use the MAAT toolbox
galCoords = celestial.coo.coco([RA dec],'j2000.0','g','r','r');
long = wrap2pi(galCoords(:,1)); 
lat = galCoords(:,2); 
end