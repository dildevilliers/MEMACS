function [az,alt] = GalLongLat2Horiz(long,lat,julDate,location)
% function [az,alt] = GalLongLat2Horiz(long,lat,julDate,location)
% Calculates horizontal coordinates from galactic coordinates.
% Inputs
% - long: Galactic longitude (radians)
% - lat: Galactic latitude (radians)
% - julDate: The Julian date (nr of days)
% - location: earth longitude and latitude in radians [long, lat]
% Outputs
% - az: Horizontal azimuth (radians)
% - alt: Horizontal altitude/elevation (radians)

[RA,dec] = GalLongLat2RAdec(long,lat);
horzCoords= wrap2pi(celestial.coo.horiz_coo([RA dec],julDate,location,'h'));
az = horzCoords(:,1);
alt = horzCoords(:,2);

end