function [az,alt] = RAdec2Horiz(RA,dec,julDate,location)
% [az,alt] = RAdec2Horiz(RA,dec,julDate,location)
% Calculates horizontal coordinates from equitorial coordinates RA and dec 
% Inputs
% - RA: right ascension (radians)
% - dec: declination (radians)
% - julDate: The Julian date (nr of days)
% - location: earth longitude and latitude in radians [long, lat]
% Outputs
% - az: Horizontal azimuth (radians)
% - alt: Horizontal altitude/elevation (radians)

% use the MAAT toolbox
horzCoords= wrap2pi(celestial.coo.horiz_coo([RA dec],julDate,location,'h'));
az = horzCoords(:,1);
alt = horzCoords(:,2);
end