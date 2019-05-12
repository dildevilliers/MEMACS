function [RA,dec] = Horiz2RAdec(az,alt,julDate,location)
% function [RA,dec] = Horiz2RAdec(az,alt,julDate,location)
% Calculates equitorial coordinates RA and dec from horizontal coordinates.
% Inputs
% - az: Horizontal azimuth (radians)
% - alt: Horizontal altitude/elevation (radians)
% - julDate: The Julian date (nr of days)
% - location: earth longitude and latitude in radians [long, lat]
% Outputs
% - RA: right ascension (radians)
% - dec: declination (radians)

% use the MAAT toolbox
equCoords= wrap2pi(celestial.coo.horiz_coo([az alt],julDate,location,'e'));
RA = equCoords(:,1); 
dec = equCoords(:,2); 
end