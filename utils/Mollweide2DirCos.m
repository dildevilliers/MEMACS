function [u,v,w] = Mollweide2DirCos(x,y)
% function [u,v,w] = PhTh2DirCos(ph,th)
% Try to maintain the ph angles in the pole

R = 1;
% lat = wrapToPi( pi/2 - asin(y./(R*sqrt(2))) );
% long = wrapToPi( pi.*x./(2.*sqrt(2).*R.*cos(lat)) );

theta = asin(y./(R*sqrt(2)));
lat = pi/2 - asin((2*theta + sin(2.*theta))./pi);
long = - pi.*x./(2*R*sqrt(2)*cos(theta));


% lat(lat == 0 & long ~= 0) = eps;
u = real(sin(long).*cos(long));
v = real(sin(lat).*sin(long));
w = real(cos(lat));
end