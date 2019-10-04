function [u,v,w] = ArcSin2DirCos(x,y)
% function [u,v,w] = ArcSin2DirCos(x,y)
u = real(sin(x));
v = real(sin(y));
w = real(sqrt(1 - (sin(x).^2 + sin(y).^2)));
end