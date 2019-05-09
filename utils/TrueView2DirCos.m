function [u,v,w] = TrueView2DirCos(x,y)
Th = sqrt(x.^2 + y.^2);
Ph = atan2(y,x);
% function [u,v,w] = TrueView2DirCos(x,y)
% Try to maintain the Ph angles in the pole
Th(Th == 0 & Ph ~= 0) = eps;
u = real(sin(Th).*cos(Ph));
v = real(sin(Th).*sin(Ph));
w = real(cos(Th));
end