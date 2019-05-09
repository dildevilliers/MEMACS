function [hHat,vHat] = unitVectorsDirCos(th,ph)
% Returns the 2 direction cosine coordinate system unit vectors as a function of
% input angles th and ph. Output vectors in Pnt3D class

hHatx = cos(th).*cos(ph).^2 + sin(ph).^2;
hHaty = 0.5.*(cos(th) - 1).*sin(2.*ph);
hHatz = -sin(th).*cos(ph);
hHat = Pnt3D(hHatx,hHaty,hHatz);

vHatx = hHaty;
vHaty = cos(th).*sin(ph).^2 + cos(ph).^2;
vHatz = -sin(th).*sin(ph);
vHat = Pnt3D(vHatx,vHaty,vHatz);


